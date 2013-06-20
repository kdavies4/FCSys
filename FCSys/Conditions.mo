within FCSys;
package Conditions "Models to specify and measure operating conditions"
  extends Modelica.Icons.Package;
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model FaceCondition "Test the conditions for the face of a subregion"
      extends Modelica.Icons.Example;

      ByConnector.FaceBus.Single.FaceBus face(gas(inclH2O=true, H2O(redeclare
              ByConnector.Face.Single.Material.Current material)))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));
      Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclFacesX=false,
        inclFacesY=true,
        inclFacesZ=false,
        inclTransX=false,
        inclTransY=true,
        graphite('inclC+'=true, 'C+'(V_IC=0.5*U.cc)),
        gas(inclH2O=true))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    equation
      connect(subregion.yPositive, face.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (experiment(NumberOfIntervals=5000), Commands(file=
              "Resources/Scripts/Dymola/Conditions.Examples.FaceCondition.mos"));
    end FaceCondition;

    model FaceConditionPhases
      "<html>Test the conditions for the face of a subregion with phases</html>"
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;
      extends Modelica.Icons.Example;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) =
        ones(3)*U.cm "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Volume V=product(L) "Volume";

      // Included components of translational momentum
      parameter Boolean inclTransX=false "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=true "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));

      // Included faces
      parameter Boolean inclFacesX=false "X" annotation (
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
      parameter Boolean inclFacesZ=false "Z" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));

      ByConnector.FaceBus.Single.Phases.Gas face(inclH2O=true, H2O(redeclare
            Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal,
            redeclare Conditions.ByConnector.Face.Single.Material.Current
            material(source(y=U.A))))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));

      Subregions.Volume volume(n_phases=1)
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      Subregions.Phases.Gas gas(
        inclH2=false,
        inclH2O=true,
        final n_faces=n_faces,
        T_IC=environment.T,
        reduceTemp=false,
        H2O(T_IC=300*U.K))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    protected
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
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
      final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
        "Translational-momentum-component indices of the Cartesian axes";
      final inner parameter Integer facesCart[Axis]=enumerate(inclFaces)
        "Face-pair indices of the Cartesian axes";

    equation
      connect(gas.yPositive, face.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(volume.inertDalton[1], gas.inertDalton) annotation (Line(
          points={{11,-11},{8,-8}},
          color={47,107,251},
          smooth=Smooth.None));
      annotation (experiment);
    end FaceConditionPhases;

    model Router
      "<html>Test the <a href=\"modelica://FCSys.Conditions.Router\">Router</a> model</html>"
      extends Modelica.Icons.Example;

      // TODO:  Make this into a meaningful example.
      Conditions.Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    end Router;

    model AnodeAdapter
      "<html>Test the <a href=\"modelica://FCSys.Conditions.Adapters.Anode\">Anode</a> adapter</html>"
      extends Modelica.Icons.Example;

      inner Modelica.Fluid.System system(T_ambient=293.15 + 5)
        annotation (Placement(transformation(extent={{40,70},{60,90}})));
      inner Conditions.Environment environment(T=350*U.K)
        annotation (Placement(transformation(extent={{70,72},{90,92}})));
      Subregions.SubregionNoIonomer subregion(
        L={1,1,1}*U.cm,
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
              "Resources/Scripts/Dymola/Conditions.Examples.Adapteminus.mos"));
    end AnodeAdapter;

  end Examples;

  package Adapters
    "<html>Adapters to the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
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

      Connectors.FaceBus face
        "Multi-species connector for translational momentum and heat"
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
            iconTransformation(extent={{-90,-10},{-70,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
          Medium = GasMedium) "Modelica fluid port for the gas" annotation (
          Placement(transformation(extent={{70,70},{90,90}}),
            iconTransformation(extent={{70,50},{90,70}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {70,30},{90,50}}), iconTransformation(extent={{70,10},{90,30}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
        "Modelica heat port" annotation (Placement(transformation(extent={{70,-14},
                {90,6}}), iconTransformation(extent={{70,-30},{90,-10}})));
      Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final package
          Medium = LiquidMedium) "Modelica fluid port for the liquid"
        annotation (Placement(transformation(extent={{70,-46},{90,-26}}),
            iconTransformation(extent={{70,-70},{90,-50}})));
      Phases.AnodeGas gas(redeclare final package Medium = GasMedium)
        "Gas subadapter"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Phases.Graphite graphite "Graphite subadapter"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
        "Liquid subadapter"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
        "Modelica translational flanges" annotation (Placement(transformation(
              extent={{70,-90},{90,-70}}), iconTransformation(extent={{70,-110},
                {90,-90}})));

    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{-8,40},{-40,40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gasPort, gas.fluidPort) annotation (Line(
          points={{80,80},{40,80},{40,44},{8,44}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(gas.heatPort, heatPort) annotation (Line(
          points={{8,36},{30,36},{30,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));

      connect(graphite.face, face.graphite) annotation (Line(
          points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));

      connect(graphite.pin, pin) annotation (Line(
          points={{8,4},{50,4},{50,40},{80,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(graphite.heatPort, heatPort) annotation (Line(
          points={{8,-4},{80,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(liquid.face, face.liquid) annotation (Line(
          points={{-8,-40},{-40,-40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquidPort, liquid.fluidPort) annotation (Line(
          points={{80,-36},{8,-36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquid.heatPort, heatPort) annotation (Line(
          points={{8,-44},{30,-44},{30,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(gas.flange, flange) annotation (Line(
          points={{8,40},{40,40},{40,-80},{80,-80}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(graphite.flange, flange) annotation (Line(
          points={{8,6.10623e-16},{24,6.10623e-16},{24,0},{40,0},{40,-80},{80,-80}},

          color={0,127,0},
          smooth=Smooth.None));

      connect(liquid.flange, flange) annotation (Line(
          points={{8,-40},{40,-40},{40,-80},{80,-80}},
          color={0,127,0},
          smooth=Smooth.None));
      annotation (Icon(graphics={Line(
                  points={{0,60},{0,-100}},
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
                  smooth=Smooth.None),Line(
                  points={{0,-100},{70,-100}},
                  color={0,127,0},
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

      Connectors.FaceBus face
        "Multi-species connector for translational momentum and heat"
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
            iconTransformation(extent={{-90,-10},{-70,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
          Medium = GasMedium) "Modelica fluid port for the gas" annotation (
          Placement(transformation(extent={{70,70},{90,90}}),
            iconTransformation(extent={{70,50},{90,70}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {70,30},{90,50}}), iconTransformation(extent={{70,10},{90,30}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
        "Modelica heat port" annotation (Placement(transformation(extent={{70,-14},
                {90,6}}), iconTransformation(extent={{70,-30},{90,-10}})));
      Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final package
          Medium = LiquidMedium) "Modelica fluid port for the liquid"
        annotation (Placement(transformation(extent={{70,-46},{90,-26}}),
            iconTransformation(extent={{70,-70},{90,-50}})));

      Phases.CathodeGas gas(redeclare final package Medium = GasMedium)
        "Gas subadapter"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Phases.Graphite graphite "Graphite subadapter"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
        "Liquid subadapter"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
        "Modelica translational flanges" annotation (Placement(transformation(
              extent={{70,-90},{90,-70}}), iconTransformation(extent={{70,-110},
                {90,-90}})));

    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{-8,40},{-40,40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gasPort, gas.fluidPort) annotation (Line(
          points={{80,80},{40,80},{40,44},{8,44}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(gas.heatPort, heatPort) annotation (Line(
          points={{8,36},{30,36},{30,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(graphite.face, face.graphite) annotation (Line(
          points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));

      connect(graphite.pin, pin) annotation (Line(
          points={{8,4},{50,4},{50,40},{80,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(graphite.heatPort, heatPort) annotation (Line(
          points={{8,-4},{80,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(liquid.face, face.liquid) annotation (Line(
          points={{-8,-40},{-40,-40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquidPort, liquid.fluidPort) annotation (Line(
          points={{80,-36},{8,-36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquid.heatPort, heatPort) annotation (Line(
          points={{8,-44},{30,-44},{30,-4},{80,-4}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(gas.flange, flange) annotation (Line(
          points={{8,40},{40,40},{40,-80},{80,-80}},
          color={0,127,0},
          smooth=Smooth.None));
      connect(graphite.flange, flange) annotation (Line(
          points={{8,6.10623e-16},{24,6.10623e-16},{24,0},{40,0},{40,-80},{80,-80}},

          color={0,127,0},
          smooth=Smooth.None));

      connect(liquid.flange, flange) annotation (Line(
          points={{8,-40},{40,-40},{40,-80},{80,-80}},
          color={0,127,0},
          smooth=Smooth.None));
      annotation (Icon(graphics={Line(
                  points={{0,60},{0,-100}},
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
                  smooth=Smooth.None),Line(
                  points={{0,-100},{70,-100}},
                  color={0,127,0},
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

        Conditions.Adapters.Species.FluidNeutral H2(redeclare package Medium =
              Modelica.Media.IdealGases.SingleGases.H2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false), redeclare package Data =
              Characteristics.H2.Gas)
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Conditions.Adapters.Species.FluidNeutral H2O(redeclare package Data =
              Characteristics.H2O.Gas (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false), redeclare final package
            Medium = Modelica.Media.IdealGases.SingleGases.H2O)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
        Junctions.Junction2 junction
          annotation (Placement(transformation(extent={{62,30},{42,50}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,30},{90,50}}), iconTransformation(
                extent={{70,30},{90,50}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));

      equation
        // H2
        connect(H2.face, face.H2) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2.heatPort, heatPort) annotation (Line(
            points={{8,16},{60,16},{60,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(H2.fluidPort, junction.purePort1) annotation (Line(
            points={{8,24},{28,24},{28,44},{44,44}},
            color={0,127,255},
            smooth=Smooth.None));

        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,-24},{60,-24},{60,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(H2O.fluidPort, junction.purePort2) annotation (Line(
            points={{8,-16},{32,-16},{32,36},{44,36}},
            color={0,127,255},
            smooth=Smooth.None));

        // Mixture
        connect(junction.mixturePort, fluidPort) annotation (Line(
            points={{60,40},{80,40}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(H2.flange, flange) annotation (Line(
            points={{8,20},{50,20},{50,5.55112e-16},{80,5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));
        connect(H2O.flange, flange) annotation (Line(
            points={{8,-20},{50,-20},{50,5.55112e-16},{80,5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,40},{70,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
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
          annotation (Placement(transformation(extent={{60,30},{40,50}})));

        Conditions.Adapters.Species.FluidNeutral H2O(redeclare package Data =
              Characteristics.H2O.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Conditions.Adapters.Species.FluidNeutral N2(redeclare package Data =
              Characteristics.N2.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.N2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Conditions.Adapters.Species.FluidNeutral O2(redeclare package Data =
              Characteristics.O2.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.O2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,30},{90,50}}), iconTransformation(
                extent={{70,30},{90,50}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));

      equation
        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.fluidPort, junction.purePort1) annotation (Line(
            points={{8,24},{26,24},{26,44},{42,44}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,16},{60,16},{60,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        // N2
        connect(N2.face, face.N2) annotation (Line(
            points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(N2.fluidPort, junction.purePort2) annotation (Line(
            points={{8,4},{30,4},{30,40},{42,40}},
            color={0,127,255},
            smooth=Smooth.None));

        connect(N2.heatPort, heatPort) annotation (Line(
            points={{8,-4},{60,-4},{60,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        // O2
        connect(O2.face, face.O2) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(O2.fluidPort, junction.purePort3) annotation (Line(
            points={{8,-16},{34,-16},{34,36},{42,36}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(O2.heatPort, heatPort) annotation (Line(
            points={{8,-24},{60,-24},{60,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        // Mixture
        connect(junction.mixturePort, fluidPort) annotation (Line(
            points={{58,40},{80,40}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(H2O.flange, flange) annotation (Line(
            points={{8,20},{50,20},{50,5.55112e-16},{80,5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));
        connect(N2.flange, flange) annotation (Line(
            points={{8,6.10623e-16},{44,6.10623e-16},{44,5.55112e-16},{80,
                5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));

        connect(O2.flange, flange) annotation (Line(
            points={{8,-20},{50,-20},{50,5.55112e-16},{80,5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,40},{70,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None)}));
      end CathodeGas;

      model Graphite
        "<html>Adapter for graphite between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        extends BaseClasses.PartialPhase;

        Species.'e-' 'e-'(redeclare package Data =
              Characteristics.'e-'.Graphite)
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));

        Species.Solid 'C+'(redeclare package Data =
              Characteristics.'C+'.Graphite)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));

      equation
        // C
        connect('C+'.face, face.'C+') annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect('C+'.heatPort, heatPort) annotation (Line(
            points={{8,-24},{40,-24},{40,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        // e-
        connect('e-'.face, face.'e-') annotation (Line(
            points={{-8,40},{-40,40},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect('e-'.heatPort, heatPort) annotation (Line(
            points={{8,36},{40,36},{40,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        connect('e-'.pin, pin) annotation (Line(
            points={{8,40},{80,40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(flange, 'C+'.flange) annotation (Line(
            points={{80,5.55112e-16},{20,5.55112e-16},{20,-20},{8,-20}},
            color={0,127,0},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,40},{70,40}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
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

        Conditions.Adapters.Species.FluidNeutral H2O(redeclare package Data =
              Characteristics.H2O.Liquid, redeclare final package Medium =
              Medium)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,30},{90,50}}), iconTransformation(
                extent={{70,30},{90,50}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));

      equation
        // H2O
        connect(H2O.face, face.H2) annotation (Line(
            points={{-8,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-80,
                5.55112e-16}},
            color={0,0,0},
            smooth=Smooth.None));

        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,-4},{40,-4},{40,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));

        connect(H2O.fluidPort, fluidPort) annotation (Line(
            points={{8,4},{40,4},{40,40},{80,40}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(flange, H2O.flange) annotation (Line(
            points={{80,5.55112e-16},{44,5.55112e-16},{44,6.10623e-16},{8,
                6.10623e-16}},
            color={0,127,0},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,40},{70,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None)}));
      end Liquid;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialPhase
          "<html>Partial adapter for a phase between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.BaseClasses.Icons.Names.Top3;

          Connectors.FaceBus face "FCSys face connector" annotation (Placement(
                transformation(extent={{-90,-10},{-70,10}}), iconTransformation(
                  extent={{-90,-10},{-70,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-50},{90,-30}}), iconTransformation(extent={{70,-50},{90,
                    -30}})));
          annotation (Icon(graphics={Line(
                          points={{0,0},{-70,0}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,-40},{70,-40}},
                          color={191,0,0},
                          smooth=Smooth.None)}));

        end PartialPhase;

      end BaseClasses;

    end Phases;

    package Species "Adapters for single species"
      extends Modelica.Icons.Package;

      model 'e-'
        "<html>Adapter to connect e<sup>-</sup> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (electrical and heat only)</html>"
        import FCSys.BaseClasses.Utilities.inSign;
        extends FCSys.BaseClasses.Icons.Names.Top2;

        // Geometry
        parameter Q.Area A=U.cm^2 "Area of the interface"
          annotation (Dialog(group="Geometry"));
        parameter Side side=Side.n "FCSys side of the interface"
          annotation (Dialog(group="Geometry"));

        replaceable package Data = Characteristics.'e-'.Graphite constrainedby
          Characteristics.BaseClasses.Characteristic
          "Characteristic data (for FCSys)" annotation (
          Dialog(group="Material properties"),
          __Dymola_choicesAllMatching=true,
          Placement(transformation(extent={{-60,40},{-40,60}}),
              iconTransformation(extent={{-10,90},{10,110}})));

        Connectors.Face face
          "Connector for material, momentum, and energy of a single species"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},{90,
                  10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -50},{90,-30}}), iconTransformation(extent={{70,-50},{90,-30}})));

      equation
        // Equal properties
        face.T = heatPort.T*U.K "Temperature";

        // Conservation (without storage)
        0 = face.Ndot "Material diffusion";
        0 = A*face.rho*face.phi[1] + pin.i*U.A/Data.z
          "Material advection (also charge)";
        inSign(side)*pin.v*face.rho*A*Data.z*U.V = face.mPhidot[1]
          "Normal translational momentum";
        {0,0} = face.mPhidot[2:3] "Transverse translational momentum";
        0 = face.Qdot + heatPort.Q_flow*U.W "Energy";
        // Note:  All of the advective terms (for all the balance equations)
        // cancel across the interface.
        annotation (Icon(graphics={Line(
                      points={{0,0},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash),Line(
                      points={{0,0},{70,0}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{70,-40}},
                      color={140,0,0},
                      smooth=Smooth.None),Line(
                      points={{-70,0},{0,0}},
                      color={127,127,127},
                      smooth=Smooth.None)}));
      end 'e-';

      model Fluid
        "<html>Adapter to connect a single fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        import FCSys.BaseClasses.Utilities.inSign;
        extends FCSys.BaseClasses.Icons.Names.Top3;

        parameter Q.Area A=U.cm^2 "Area of the interface"
          annotation (Dialog(group="Geometry"));
        parameter Side side=Side.n "FCSys side of the interface"
          annotation (Dialog(group="Geometry"));
        replaceable package Data = Characteristics.BaseClasses.Characteristic
          "Characteristic data (for FCSys)" annotation (
          Dialog(group="Material properties"),
          __Dymola_choicesAllMatching=true,
          Placement(transformation(extent={{-60,40},{-40,60}}),
              iconTransformation(extent={{-10,90},{10,110}})));
        replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2O
          constrainedby Modelica.Media.Interfaces.PartialPureSubstance
          "Medium model (for Modelica)" annotation (choicesAllMatching=true,
            Dialog(group="Material properties"));

        Medium.BaseProperties medium "Base properties of the fluid";
        Q.Current I "Material current";

        Connectors.Face face
          "Connector for material, momentum, and energy of a single species"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,30},{90,50}}), iconTransformation(
                extent={{70,30},{90,50}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -50},{90,-30}}), iconTransformation(extent={{70,-50},{90,-30}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,-90},{90,-70}}), iconTransformation(extent={{70,-90},{90,
                  -70}})));

      equation
        // Aliases (for common terms)
        I = face.Ndot + face.phi[1]*face.rho*A "Current";

        // Base media properties
        medium.p = fluidPort.p;
        medium.T = heatPort.T;
        medium.Xi = ones(Medium.nXi)/Medium.nXi;

        // Equal properties
        medium.MM*face.rho = medium.d*U.mol/U.m^3 "Density";
        face.phi = der(flange.s)*U.m/U.s "Velocity";
        face.T = heatPort.T*U.K "Temperature";
        medium.h = fluidPort.h_outflow;

        // Conservation (without storage)
        0 = Data.z*I + pin.i*U.A "Charge";
        0 = I + (fluidPort.m_flow/medium.MM)*U.mol/U.s "Material";
        inSign(side)*pin.v*face.rho*A*Data.z*U.V = face.mPhidot[1] + flange[1].f
          *U.N "Normal translational momentum";
        {0,0} = face.mPhidot[2:3] + flange[2:3].f*U.N
          "Transverse translational momentum";
        0 = face.Qdot + heatPort.Q_flow*U.W "Energy";
        // Note:  All of the advective terms (for all the balance equations)
        // cancel across the interface.
        annotation (Icon(graphics={Line(
                      points={{0,40},{70,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,40},{0,-80}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None),Line(
                      points={{0,-80},{70,-80}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{70,-40}},
                      color={140,0,0},
                      smooth=Smooth.None),Line(
                      points={{-70,0},{0,0}},
                      color={127,127,127},
                      smooth=Smooth.None)}));
      end Fluid;

      model FluidNeutral
        "<html>Adapter to connect a single neutral fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        import assert = FCSys.BaseClasses.Utilities.assertEval;
        extends FCSys.BaseClasses.Icons.Names.Top3;

        parameter Q.Area A=U.cm^2 "Area of the interface"
          annotation (Dialog(group="Geometry"));
        replaceable package Data = Characteristics.BaseClasses.Characteristic
          "Characteristic data (for FCSys)" annotation (
          Dialog(group="Material properties"),
          __Dymola_choicesAllMatching=true,
          Placement(transformation(extent={{-60,40},{-40,60}}),
              iconTransformation(extent={{-10,90},{10,110}})));
        replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2O
          constrainedby Modelica.Media.Interfaces.PartialPureSubstance
          "Medium model (for Modelica)" annotation (choicesAllMatching=true,
            Dialog(group="Material properties"));

        Medium.BaseProperties medium "Base properties of the fluid";

        Connectors.Face face
          "Connector for material, momentum, and energy of a single species"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,30},{90,50}}), iconTransformation(
                extent={{70,30},{90,50}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}),iconTransformation(extent={{70,-10},
                  {90,10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -50},{90,-30}}), iconTransformation(extent={{70,-50},{90,-30}})));

      initial equation
        assert(Data.z == 0,
          "The species must be neutral, but its chemical formula is " + Data.formula);

      equation
        // Base media properties
        medium.p = fluidPort.p;
        medium.T = heatPort.T;
        medium.Xi = ones(Medium.nXi)/Medium.nXi;

        // Equal properties
        medium.MM*face.rho = medium.d*U.mol/U.m^3 "Density";
        face.phi = der(flange.s)*U.m/U.s "Velocity";
        face.T = heatPort.T*U.K "Temperature";
        medium.h = fluidPort.h_outflow;

        // Conservation (without storage)
        0 = face.Ndot + face.phi[1]*face.rho*A + (fluidPort.m_flow/medium.MM)*U.mol
          /U.s "Material";
        {0,0,0} = face.mPhidot + flange.f*U.N "Translational momentum";
        0 = face.Qdot + heatPort.Q_flow*U.W "Energy";
        // Note:  All of the advective terms (for all the balance equations)
        // cancel across the interface.
        annotation (Icon(graphics={Line(
                      points={{0,40},{70,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{70,-40}},
                      color={140,0,0},
                      smooth=Smooth.None),Line(
                      points={{-70,0},{0,0}},
                      color={127,127,127},
                      smooth=Smooth.None)}));
      end FluidNeutral;

      model Solid
        "<html>Adapter to connect a single solid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (heat only)</html>"
        extends FCSys.BaseClasses.Icons.Names.Top2;

        replaceable package Data = Characteristics.BaseClasses.Characteristic
          "Characteristic data (for FCSys)" annotation (
          Dialog(group="Material properties"),
          __Dymola_choicesAllMatching=true,
          Placement(transformation(extent={{-60,40},{-40,60}}),
              iconTransformation(extent={{-10,90},{10,110}})));

        Connectors.Face face
          "Connector for material, momentum, and energy of a single species"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges"
          annotation (Placement(transformation(extent={{70,-10},{90,10}})));

        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -50},{90,-30}}), iconTransformation(extent={{70,-50},{90,-30}})));

      equation
        // Equal properties
        face.phi = der(flange.s)*U.m/U.s "Velocity";
        face.T = heatPort.T*U.K "Temperatures";

        // Conservation (without storage)
        0 = face.Ndot "Material";
        {0,0,0} = face.mPhidot + flange.f*U.N "Translational momentum";
        0 = face.Qdot + heatPort.Q_flow*U.W "Energy";
        // Note:  All of the advective terms (for all the balance equations)
        // cancel across the interface.
        annotation (Icon(graphics={Line(
                      points={{0,0},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None),Line(
                      points={{-70,0},{0,0}},
                      color={127,127,127},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{70,-40}},
                      color={140,0,0},
                      smooth=Smooth.None)}));
      end Solid;

    end Species;

    package Junctions
      "<html><a href=\"modelica://Modelica\">Modelica</a> junctions between pure substances and their mixtures</html>"
      extends Modelica.Icons.Package;

      model Junction2 "Junction between two pure substances and their mixture"
        import assert = FCSys.BaseClasses.Utilities.assertEval;
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
          annotation (Placement(transformation(extent={{70,30},{90,50}}),
              iconTransformation(extent={{70,30},{90,50}})));
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
          Documentation(info="<html><p>Assumptions:
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol></p></html>"),
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
        import assert = FCSys.BaseClasses.Utilities.assertEval;
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
          annotation (Placement(transformation(extent={{70,30},{90,50}}),
              iconTransformation(extent={{70,30},{90,50}})));
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
          Documentation(info="<html><p>Assumptions:
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol></p></html>"),
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

      package BaseClasses "Base classes (generally not for direct use)"
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

          Modelica.SIunits.MassFraction X[MixtureMedium.nX]
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
      "<html><a href=\"modelica://Modelica.Media\">Modelica media</a> models to interface with the <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">cell</a></html>"
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
      extends Modelica.Icons.UnderConstruction;
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
      Connectors.RealOutputBus y "Output signals as a bus" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-110}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-50})));

    protected
      Modelica.Blocks.Math.Add sumAnMFC annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,30})));
      Modelica.Blocks.Math.Add sumCaMFC annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,30})));

      Modelica.Blocks.Math.Gain from_V(k=U.V) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-152,0})));
      Modelica.Blocks.Math.Gain 'from1_%'(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-134,-30})));
      Modelica.Blocks.Math.Gain 'from2_%'(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-114,0})));
      Blocks.UnitConversions.From_kPag from1_kPag annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-94,-30})));
      Blocks.UnitConversions.From_kPag from2_kPag annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-74,0})));
      Blocks.UnitConversions.From_kPag from3_kPag annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-54,-30})));
      Blocks.UnitConversions.From_kPag from4_kPag annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-34,0})));
      Modelica.Blocks.Math.Gain from1_LPM(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,-30})));
      Modelica.Blocks.Math.Gain from2_LPM(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,0})));
      Blocks.UnitConversions.From_degC from1_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={32,-30})));
      Blocks.UnitConversions.From_degC from2_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,0})));
      Blocks.UnitConversions.From_degC from3_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={72,-30})));
      Blocks.UnitConversions.From_degC from4_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={92,0})));
      Blocks.UnitConversions.From_degC from5_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={112,-30})));
      Blocks.UnitConversions.From_degC from6_degC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={132,0})));
      Modelica.Blocks.Math.Gain from_A(k=U.A) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={152,-30})));

      Connectors.RealOutputInternal v(final unit="l2.m/(N.T2)")
        "CVM Cell 1 Voltage" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-152,-60})));
      Connectors.RealOutputInternal RHAnFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Anode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-134,-60})));
      Connectors.RealOutputInternal RHCaFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Cathode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-114,-60})));
      Connectors.RealOutputInternal p_anFPNegY(final unit="m/(l.T2)")
        "Pressure anode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-94,-60})));
      Connectors.RealOutputInternal p_anFPPosY(final unit="m/(l.T2)", final min
          =0) "Pressure anode outlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-74,-60})));
      Connectors.RealOutputInternal p_caFPNegY(final unit="m/(l.T2)")
        "Pressure cathode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-54,-60})));
      Connectors.RealOutputInternal p_caFPPosY(final unit="m/(l.T2)", final min
          =0) "Pressure anode outlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-34,-60})));
      Connectors.RealOutputInternal Vdot_anFPNegY_H2(final unit="l3/T")
        "Flow anode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,-60})));
      Connectors.RealOutputInternal Vdot_caFPNegY_air(final unit="l3/T")
        "Flow cathode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,-60})));
      Connectors.RealOutputInternal T_anFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={32,-60})));
      Connectors.RealOutputInternal T_anFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,-60})));
      Connectors.RealOutputInternal T_caFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={72,-60})));
      Connectors.RealOutputInternal T_caFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={92,-60})));
      Connectors.RealOutputInternal T_anFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate anode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={112,-60})));
      Connectors.RealOutputInternal T_caFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate cathode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={132,-60})));
      Connectors.RealOutputInternal I(final unit="N/T") "Measured load"
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
      connect(from_V.u, combiTimeTable.y[1]) annotation (Line(
          points={{-152,12},{-152,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect('from1_%'.u, combiTimeTable.y[2]) annotation (Line(
          points={{-134,-18},{-134,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect('from2_%'.u, combiTimeTable.y[3]) annotation (Line(
          points={{-114,12},{-114,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(from1_kPag.u, combiTimeTable.y[4]) annotation (Line(
          points={{-94,-19},{-94,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from2_kPag.u, combiTimeTable.y[5]) annotation (Line(
          points={{-74,11},{-74,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from3_kPag.u, combiTimeTable.y[6]) annotation (Line(
          points={{-54,-19},{-54,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from4_kPag.u, combiTimeTable.y[7]) annotation (Line(
          points={{-34,11},{-34,50},{-1.44329e-15,50},{-1.44329e-15,69}},
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
      connect(from1_degC.u, combiTimeTable.y[12]) annotation (Line(
          points={{32,-19},{32,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from2_degC.u, combiTimeTable.y[13]) annotation (Line(
          points={{52,11},{52,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from3_degC.u, combiTimeTable.y[14]) annotation (Line(
          points={{72,-19},{72,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from4_degC.u, combiTimeTable.y[15]) annotation (Line(
          points={{92,11},{92,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from5_degC.u, combiTimeTable.y[16]) annotation (Line(
          points={{112,-19},{112,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from6_degC.u, combiTimeTable.y[17]) annotation (Line(
          points={{132,11},{132,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(from_A.u, combiTimeTable.y[18]) annotation (Line(
          points={{152,-18},{152,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from unit conversion to internal outputs
      connect(v, from_V.y) annotation (Line(
          points={{-152,-60},{-152,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHAnFPNegX, 'from1_%'.y) annotation (Line(
          points={{-134,-60},{-134,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHCaFPNegX, 'from2_%'.y) annotation (Line(
          points={{-114,-60},{-114,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPNegY, from1_kPag.y) annotation (Line(
          points={{-94,-60},{-94,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPPosY, from2_kPag.y) annotation (Line(
          points={{-74,-60},{-74,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPNegY, from3_kPag.y) annotation (Line(
          points={{-54,-60},{-54,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPPosY, from4_kPag.y) annotation (Line(
          points={{-34,-60},{-34,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_anFPNegY_H2, from1_LPM.y) annotation (Line(
          points={{-14,-60},{-14,-50.5},{-14,-41},{-14,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_caFPNegY_air, from2_LPM.y) annotation (Line(
          points={{14,-60},{14,-35.5},{14,-11},{14,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, from1_degC.y) annotation (Line(
          points={{32,-60},{32,-50.5},{32,-41},{32,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, from2_degC.y) annotation (Line(
          points={{52,-60},{52,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, from3_degC.y) annotation (Line(
          points={{72,-60},{72,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, from4_degC.y) annotation (Line(
          points={{92,-60},{92,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, from5_degC.y) annotation (Line(
          points={{112,-60},{112,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPX, from6_degC.y) annotation (Line(
          points={{132,-60},{132,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(I, from_A.y) annotation (Line(
          points={{152,-60},{152,-41}},
          color={0,0,127},
          smooth=Smooth.None));

      // Summations
      connect(sumAnMFC.y, from1_LPM.u) annotation (Line(
          points={{-14,19},{-14,9.75},{-14,9.75},{-14,0.5},{-14,-18},{-14,-18}},

          color={0,0,127},
          smooth=Smooth.None));

      connect(sumCaMFC.y, from2_LPM.u) annotation (Line(
          points={{14,19},{14,17.25},{14,17.25},{14,15.5},{14,12},{14,12}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from internal outputs to public output
      connect(v, y.v) annotation (Line(
          points={{-152,-60},{-152,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHAnFPNegX, y.RHAnFPNegX) annotation (Line(
          points={{-134,-60},{-134,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHCaFPNegX, y.RHCaFPNegX) annotation (Line(
          points={{-114,-60},{-114,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPNegY, y.p_anFPNegY) annotation (Line(
          points={{-94,-60},{-94,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPPosY, y.p_anFPPosY) annotation (Line(
          points={{-74,-60},{-74,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPNegY, y.p_caFPNegY) annotation (Line(
          points={{-54,-60},{-54,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPPosY, y.p_caFPPosY) annotation (Line(
          points={{-34,-60},{-34,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_anFPNegY_H2, y.Vdot_anFPNegY_H2) annotation (Line(
          points={{-14,-60},{-14,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_caFPNegY_air, y.Vdot_caFPNegY_air) annotation (Line(
          points={{14,-60},{14,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, y.T_anFPNegY) annotation (Line(
          points={{32,-60},{32,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, y.T_anFPPosY) annotation (Line(
          points={{52,-60},{52,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, y.T_caFPNegY) annotation (Line(
          points={{72,-60},{72,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, y.T_caFPPosY) annotation (Line(
          points={{92,-60},{92,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, y.T_anFPX) annotation (Line(
          points={{112,-60},{112,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T_caFPX, y.T_caFPX) annotation (Line(
          points={{132,-60},{132,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(I, y.I) annotation (Line(
          points={{152,-60},{152,-80},{5.55112e-16,-80},{5.55112e-16,-110}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-180,
                -100},{180,100}}), graphics), experiment(StopTime=15481,
            Algorithm="Euler"));
    end Replay;

    package BaseClasses "Base classes (generally not for direct use)"
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
          annotation (HideResult=true, choices(__Dymola_checkBox=true));

        Connectors.FaceBus anEnd[n_y, n_z] "Anode end plate" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-160,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-160,0})));
        Connectors.FaceBus caEnd[n_y, n_z] "Cathode end plate" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={160,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={162,0})));
        Connectors.FaceBus anSource[n_x_an, n_z] "Anode source" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,-160})));
        Connectors.FaceBus anSink[n_x_an, n_z] "Anode sink" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,160})));
        Connectors.FaceBus caSource[n_x_ca, n_z] "Cathode source" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-160}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,160})));
        Connectors.FaceBus caSink[n_x_ca, n_z] "Cathode sink" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,-160})));

        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anEndBC[n_y, n_z]
          (each graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-136,0})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caEndBC[n_y, n_z]
          (each graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={136,0})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anSourceBC[n_x_an,
          n_z](each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-136})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anSinkBC[n_x_an,
          n_z](each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-40,136})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caSourceBC[n_x_ca,
          n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-136})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caSinkBC[n_x_ca,
          n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,136})));
        Connectors.RealInputBus u[n_y, n_z] if inclIO annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={-160,160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={-166,166})));
        Connectors.RealOutputBus y[n_y, n_z] if inclIO annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={160,-160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={166,-166})));
        replaceable Conditions.ByConnector.FaceBus.Pair.FaceBus current[n_y,
          n_z](graphite('inclC+'=true, 'incle-'=true)) if inclIO constrainedby
          Conditions.ByConnector.FaceBus.Pair.FaceBus(graphite('inclC+'=true,
              'incle-'=true))
          annotation (Placement(transformation(extent={{-140,20},{-120,40}})));

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
            points={{-160,160},{-130,130},{-130,35}},
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
        extends Modelica.Icons.UnderConstruction;

        final parameter Integer n_x_an=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x an</sub>)</html>";
        final parameter Integer n_x_ca=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x ca</sub>)</html>";
        final parameter Integer n_y=1
          "<html>Number of subregions along the channel (<i>n</i><sub>y</sub>)</html>";
        final parameter Integer n_z=1
          "<html>Number of subregions across the channel (<i>n</i><sub>z</sub>)</html>";

        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anEnd[n_y, n_z](
            each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Normal.CurrentAreic normal(redeclare
                  Modelica.Blocks.Sources.Ramp source(height=U.A/U.cm^2,
                    duration=50))))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-30,0})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caEnd[n_y, n_z](
            each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Normal.CurrentAreic normal(redeclare
                  Modelica.Blocks.Sources.Ramp source(height=U.A/U.cm^2,
                    duration=50))))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,0})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anSource[n_x_an,
          n_z](each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-30})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated anSink[n_x_an,
          n_z](each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,30})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caSource[n_x_ca,
          n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,-30})));
        Conditions.ByConnector.FaceBus.Single.FaceBusIsolated caSink[n_x_ca,
          n_z](each gas(
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

  package ByConnector "Conditions for each type of connector"
    extends Modelica.Icons.Package;
    package ChemicalReaction
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.ChemicalReaction\">ChemicalReaction</a> connector</html>"
      extends Modelica.Icons.Package;

      model ChemicalReaction
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.ChemicalReaction\">ChemicalReaction</a> connector, with efforts by default</html>"
        import FCSys.BaseClasses.Utilities.countTrue;
        import FCSys.BaseClasses.Utilities.enumerate;
        import FCSys.BaseClasses.Utilities.index;
        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        // Conditions
        replaceable Material.ReactionRate material constrainedby
          Conditions.ByConnector.ChemicalReaction.Material.BaseClasses.PartialCondition
          "Material" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{-80,46},{-60,66}})));
        replaceable Translational.Velocity translationalX(final axis) if
          inclTransX constrainedby
          Conditions.ByConnector.ChemicalReaction.Translational.BaseClasses.PartialCondition(
            axis=Axis.x) "X-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransX),
          Placement(transformation(extent={{-52,32},{-32,52}})));
        replaceable Translational.Velocity translationalY(final axis) if
          inclTransY constrainedby
          Conditions.ByConnector.ChemicalReaction.Translational.BaseClasses.PartialCondition(
            axis=Axis.y) "Y-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransY),
          Placement(transformation(extent={{-24,18},{-4,38}})));
        replaceable Translational.Velocity translationalZ(final axis) if
          inclTransZ constrainedby
          Conditions.ByConnector.ChemicalReaction.Translational.BaseClasses.PartialCondition(
            axis=Axis.z) "Z-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransZ),
          Placement(transformation(extent={{4,4},{24,24}})));
        replaceable ThermalAdvection.SpecificEntropyTemperature
          thermalAdvection(source(y=3000*U.K)) constrainedby
          Conditions.ByConnector.ChemicalReaction.ThermalAdvection.BaseClasses.PartialCondition
          "Thermal advection" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{32,-10},{52,10}})));
        replaceable ThermalDiffusion.Temperature thermalDiffusion(source(y=
                298.15*U.K)) constrainedby
          Conditions.ByConnector.ChemicalReaction.ThermalDiffusion.BaseClasses.PartialCondition
          "Thermal diffusion" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{60,-24},{80,-4}})));

        Connectors.ChemicalReaction chemical(final n_trans=n_trans)
          "Connector for a chemical reaction"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        Connectors.RealInputBus u
          "Input bus for values of specified conditions" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealOutputBus y "Output bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer cartTrans[n_trans]=index({inclTransX,
            inclTransY,inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

      equation
        // Material
        connect(material.chemical, chemical) annotation (Line(
            points={{-70,52},{-70,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.material, material.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,56},{-81,56}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(material.y, y.material) annotation (Line(
            points={{-59,56},{90,56},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // X-axis translational
        connect(translationalX.chemical, chemical) annotation (Line(
            points={{-42,38},{-42,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalX, translationalX.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,42},{-53,42}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalX.y, y.translationalX) annotation (Line(
            points={{-31,42},{90,42},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Y-axis translational
        connect(translationalY.chemical, chemical) annotation (Line(
            points={{-14,24},{-14,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalY, translationalY.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,28},{-25,28}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalY.y, y.translationalY) annotation (Line(
            points={{-3,28},{90,28},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Z-axis translational
        connect(translationalZ.chemical, chemical) annotation (Line(
            points={{14,10},{14,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalZ, translationalZ.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,14},{3,14}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalZ.y, y.translationalZ) annotation (Line(
            points={{25,14},{90,14},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Thermal advection
        connect(thermalAdvection.chemical, chemical) annotation (Line(
            points={{42,-4},{42,-28},{0,-28},{0,-40},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.thermalAdvection, thermalAdvection.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,6.10623e-16},{31,
                6.10623e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));

        connect(thermalAdvection.y, y.thermalAdvection) annotation (Line(
            points={{53,6.10623e-16},{90,6.10623e-16},{90,5.55112e-16},{110,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Thermal diffusion
        connect(thermalDiffusion.chemical, chemical) annotation (Line(
            points={{70,-18},{70,-28},{0,-28},{0,-40},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.thermalDiffusion, thermalDiffusion.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,-14},{59,-14}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(thermalDiffusion.y, y.thermalDiffusion) annotation (Line(
            points={{81,-14},{90,-14},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));
        annotation (defaultComponentName="chemical");
      end ChemicalReaction;

      model ChemicalReactionNoFlow
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.ChemicalReaction\">ChemicalReaction</a> connector, with zero flows by default</html>"
        extends ChemicalReaction(
          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.Material.Potential material,

          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.Translational.Force
            translationalX(final axis),
          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.Translational.Force
            translationalY(final axis),
          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.Translational.Force
            translationalZ(final axis),
          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.ThermalAdvection.HeatRate
            thermalAdvection,
          redeclare replaceable
            Conditions.ByConnector.ChemicalReaction.ThermalDiffusion.HeatRate
            thermalDiffusion);

        // Note:  Dymola 7.4 requires that the redeclared models are
        // resolved to the root of the library, e.g.,
        //   "Conditions.ByConnector.ChemicalReaction.Material.Potential"
        // instead of
        //   "Material.Potential".
        // Otherwise the following error is given:
        //   "Cannot show paramter [sic] dialog for redeclared class [...] since:
        //   Could not find type of redeclare :[...] in scope [...]"
        // Similar notes apply to the other top-level extended models in
        // FCSys.Conditions.ByConnector.
        annotation (defaultComponentName="chemical");

      end ChemicalReactionNoFlow;

      package Material "Conditions for additivity of volume"
        extends Modelica.Icons.Package;

        model ReactionRate
          "Specify current (measure electrochemical potential)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.ReactionRate,
            u(final unit="N/T"),
            final y(final unit="l2.m/(N.T2)") = chemical.mu);

        equation
          chemical.Ndot = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="material");
        end ReactionRate;

        model Potential "Specify electrochemical potential (measure current)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Potential,
            u(final unit="l2.m/(N.T2)"),
            final y(final unit="N/T") = chemical.Ndot);

        equation
          chemical.mu = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="material");
        end Potential;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=chemical.Ndot);

          Real x=chemical.mu "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="material",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.mu</code> and/or <code>chemical.Ndot</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model of a material condition"

            extends ByConnector.ChemicalReaction.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          equation
            // Zero values of other flows
            chemical.mPhidot = zeros(n_trans) "Force";
            chemical.Qdot_A = 0 "Rate of thermal advection";
            chemical.Qdot_D = 0 "Rate of thermal diffusion";
            annotation (defaultComponentName="material");
          end PartialCondition;

          type ConditionType = enumeration(
              ReactionRate
                "Specify reaction rate (measure electrochemical potential)",
              Potential
                "Specify electrochemical potential (measure reaction rate)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Material;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;
        model Velocity "Specify velocity (measure force)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"),
            final y(final unit="l.m/T2") = chemical.mPhidot[transCart[axis]]);

        equation
          if n_trans > 0 then
            chemical.phi[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Velocity;

        model Force "Specify force (measure velocity)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l.m/T2"),
            final y(final unit="l/T") = chemical.phi[transCart[axis]]);

        equation
          if n_trans > 0 then
            chemical.mPhidot[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Force;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=chemical.mPhidot[transCart[axis]]);

          Real x=chemical.phi[transCart[axis]]
            "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          if n_trans > 0 then
            x = u_final;
          end if;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="translational",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.phi[transCart[axis]]</code> and/or <code>chemical.mPhidot[transCart[axis]]</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a translational condition"
            extends ByConnector.ChemicalReaction.BaseClasses.PartialCondition;

            parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          protected
            outer parameter Integer cartTrans[:]
              "Cartesian-axis indices of the components of translational momentum";
            outer parameter Integer transCart[Axis]
              "Translational-momentum-component indices of the Cartesian axes";

          equation
            // Zero values of other flows
            chemical.mu = 0 "Electrochemical potential";
            for i in 1:n_trans loop
              if cartTrans[i] <> axis then
                chemical.mPhidot[i] = 0 "Force along the other axes";
              end if;
            end for;
            chemical.Qdot_A = 0 "Rate of thermal advection";
            chemical.Qdot_D = 0 "Rate of thermal diffusion";
            annotation (defaultComponentName="translational");
          end PartialCondition;

          type ConditionType = enumeration(
              Velocity "Specify velocity (measure force)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Translational;

      package ThermalAdvection "Conditions for thermal advection"
        extends Modelica.Icons.Package;

        model SpecificEntropyTemperature
          "Specify specific entropy-temperature product (measure heat flow rate)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.SpecificEntropyTemperature,

            u(final unit="l2.m/(N.T2)"),
            final y(final unit="l2.m/T3") = chemical.Qdot_A);

        equation
          chemical.sT = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermalAdvection");
        end SpecificEntropyTemperature;

        model HeatRate
          "Specify heat flow rate (measure specific entropy-temperature product)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.HeatRate,
            u(final unit="l2.m/T3"),
            final y(final unit="l2.m/(N.T2)") = chemical.sT);

        equation
          chemical.Qdot_A = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermalAdvection");
        end HeatRate;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=chemical.Qdot_A);

          Real x=chemical.sT "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="thermalAdvection",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.sT</code> and/or <code>chemical.Qdot_A</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition "Partial model for a fluid condition"

            extends ByConnector.ChemicalReaction.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          equation
            // Zero values of other flows
            chemical.mu = 0 "Electrochemical potential";
            chemical.mPhidot = zeros(n_trans) "Force";
            chemical.Qdot_D = 0 "Rate of thermal diffusion";
            annotation (defaultComponentName="thermalAdvection");
          end PartialCondition;

          type ConditionType = enumeration(
              SpecificEntropyTemperature
                "Specify specific entropy-temperature product (measure heat flow rate)",

              HeatRate
                "Specify heat flow rate (measure specific entropy-temperature product)",

              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end ThermalAdvection;

      package ThermalDiffusion "Conditions for thermal diffusion"
        extends Modelica.Icons.Package;

        model Temperature "Specify temperature (measure heat flow rate)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"),
            final y(final unit="l2.m/T3") = chemical.Qdot_D);

        equation
          chemical.T = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermalDiffusion");
        end Temperature;

        model HeatRate "Specify heat flow rate (measure temperature)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.HeatRate,
            u(final unit="l2.m/T3"),
            final y(final unit="l2.m/(N.T2)") = chemical.T);

        equation
          chemical.Qdot_D = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermalDiffusion");
        end HeatRate;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=chemical.Qdot_D);

          Real x=chemical.T "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="thermalDiffusion",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.T</code> and/or <code>chemical.Qdot_D</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition "Partial model for a fluid condition"

            extends ByConnector.ChemicalReaction.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          equation
            // Zero values of other flows
            chemical.mu = 0 "Electrochemical potential";
            chemical.mPhidot = zeros(n_trans) "Force";
            chemical.Qdot_A = 0 "Rate of thermal advection";
            annotation (defaultComponentName="thermalDiffusion");
          end PartialCondition;

          type ConditionType = enumeration(
              Temperature "Specify temperature (measure heat flow rate)",
              HeatRate "Specify heat flow rate (measure temperature)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end ThermalDiffusion;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialCondition "Partial model of a condition"
          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification",enable=internal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));
          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));
          Connectors.ChemicalReaction chemical(final n_trans=n_trans)
            "Connector for a chemical reaction"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        protected
          outer parameter Integer n_trans
            "Number of components of translational momentum";

          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                  5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));

        end PartialCondition;

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={239,142,1},
              fillPattern=FillPattern.Solid,
              fillColor={255,195,38})}));

    end ChemicalReaction;

    package ChemicalSpecies
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.ChemicalSpecies\">ChemicalSpecies</a> connector</html>"
      extends Modelica.Icons.Package;

      model Potential "Specify chemical potential (measure current)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Potential,
          u(final unit="l2.m/(N.T2)"),
          final y(final unit="N/T") = chemical.Ndot);

      equation
        chemical.mu = u_final;

      end Potential;

      model Current "Specify current (measure chemical potential)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Current,
          u(final unit="N/T"),
          final y(final unit="l2.m/(N.T2)") = chemical.mu);

      equation
        chemical.Ndot = u_final;

      end Current;

      model Custom "Custom"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=chemical.Ndot);

        Real x=chemical.mu "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (defaultComponentName="chemical", Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.mu</code> and/or <code>chemical.Ndot</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model of a material condition"
          import FCSys.BaseClasses.Utilities.countTrue;
          import FCSys.BaseClasses.Utilities.index;
          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification of material condition"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification of material condition",enable=internal),

            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));

          // Properties upon outflow
          parameter Q.Velocity phi[Axis]={0,0,0}
            "<html>Velocity (<b>&phi;</b>)</html>"
            annotation (Dialog(group="Properties upon outflow"));
          parameter Q.PotentialAbsolute sT(start=3000*U.K)
            "Specific entropy-temperature product"
            annotation (Dialog(group="Properties upon outflow"));

          // Included components of translational momentum
          parameter Boolean inclTransX=true "X" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransY=true "Y" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransZ=true "Z" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with
          // the results.

          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));
          output Q.Velocity phi_actual[n_trans]=actualStream(chemical.phi)
            "Velocity of the actual stream";
          output Q.Potential sT_actual=actualStream(chemical.sT)
            "Specific entropy-temperature product of the actual stream";
          Connectors.ChemicalSpecies chemical(final n_trans=n_trans)
            "Connector for a species of a chemical reaction"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        protected
          final parameter Integer n_trans=countTrue({inclTransX,inclTransY,
              inclTransZ}) "Number of components of translational momentum";
          final parameter Integer cartTrans[n_trans]=index({inclTransX,
              inclTransY,inclTransZ})
            "Cartesian-axis indices of the components of translational momentum";
          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          chemical.phi = phi[cartTrans];
          chemical.sT = sT;

          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                  5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (defaultComponentName="chemical");
        end PartialCondition;

        type ConditionType = enumeration(
            Current "Specify current (measure chemical potential)",
            Potential "Specify chemical potential (measure current)",
            Custom "Custom") "Types of conditions";

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={239,142,1},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}), Ellipse(
              extent={{-40,20},{20,-40}},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None)}));

    end ChemicalSpecies;

    package PhysicalBus
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.PhysicalBus\">PhysicalBus</a> or <a href=\"modelica://FCSys.Connectors.PhysicalBusInternal\">PhysicalBusInternal</a> connector</html>"
      extends Modelica.Icons.Package;

      model PhysicalBus
        "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.PhysicalBus\">PhysicalBus</a> connector, with effort by default</html>"

        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Conditionally include species.
        parameter Boolean 'inclC+'=false
          "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Physical.Potential 'C+'(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if 'inclC+' constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="C+") "Conditions" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC+'),
          Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'inclC19HF37O5S-'=false
          "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Physical.Potential 'C19HF37O5S-'(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if 'inclC19HF37O5S-' constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="inclC19HF37O5S-") "Conditions" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC19HF37O5S-'),
          Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Physical.Potential 'e-'(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if 'incle-' constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="e-") "Conditions" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'),
          Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Physical.Potential 'H+'(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if 'inclH+' constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="H+") "Conditions" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'),
          Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Physical.Potential H2(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if inclH2 constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="H2") "Conditions" annotation (
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

        replaceable Physical.Potential H2O(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if inclH2O constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="H2O") "Conditions" annotation (
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

        replaceable Physical.Potential N2(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if inclN2 constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="N2") "Conditions" annotation (
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

        replaceable Physical.Potential O2(
          final inclTransX,
          final inclTransY,
          final inclTransZ,
          final formula) if inclO2 constrainedby
          Conditions.ByConnector.Physical.BaseClasses.PartialCondition(
          inclTransX=inclTransX,
          inclTransY=inclTransY,
          inclTransZ=inclTransZ,
          formula="O2") "Conditions" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2),
          Placement(transformation(extent={{-10,-10},{10,10}})));

        // Assumptions
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        Connectors.PhysicalBus physical "Bus of multiple species"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        Connectors.RealInputBus u "Bus of inputs to specify conditions"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));

        Connectors.RealOutputBus y "Bus of measurement outputs" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

      equation
        // Note:  It would be helpful if Modelica allowed elements of expandable
        // connectors to be named by the contents of a string variable and the
        // name of an instance of a model was accessible through a string (like
        // %name is expanded to be the name of the instance of the model).  Then,
        // the connection equations that follow could be generic.

        // C+
        connect('C+'.physical, physical.'C+') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.'C+', 'C+'.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C+'.y, y.'C+') annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S-
        connect('C19HF37O5S-'.physical, physical.'C19HF37O5S-') annotation (
            Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.y, y.'C19HF37O5S-') annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.physical, physical.'e-') annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('e-'.y, y.'e-') annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.physical, physical.'H+') annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('H+'.y, y.'H+') annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2
        connect(H2.physical, physical.H2) annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.H2, H2.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2.y, y.H2) annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.physical, physical.H2O) annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.physical, physical.N2) annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.N2, N2.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(N2.y, y.N2) annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.physical, physical.O2) annotation (Line(
            points={{6.10623e-16,-4},{5.55112e-16,-40}},
            color={38,196,52},
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-11,0},{-11,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(O2.y, y.O2) annotation (Line(
            points={{11,6.10623e-16},{11,0},{110,0},{110,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (defaultComponentName="physical");
      end PhysicalBus;

      model PhysicalBusIsolated
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.PhysicalBus\">PhysicalBus</a> connector, with zero flows by default</html>"
        extends PhysicalBus(
          redeclare replaceable Conditions.ByConnector.Physical.Current 'C+'(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current
            'C19HF37O5S-'(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current 'e-'(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current 'H+'(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current H2(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current H2O(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current N2(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula),
          redeclare replaceable Conditions.ByConnector.Physical.Current O2(
            final inclTransX,
            final inclTransY,
            final inclTransZ,
            final formula));
        annotation (defaultComponentName="physical");

      end PhysicalBusIsolated;
      annotation (Documentation(info="<html><p>All of the submodels for the individual species in
  <a href=\"modelica://FCSys.Conditions.ByConnector.ChemicalBus.Gas\">Gas</a>,
  <a href=\"modelica://FCSys.Conditions.ByConnector.ChemicalBus.Graphite\">Graphite</a>,
  <a href=\"modelica://FCSys.Conditions.ByConnector.ChemicalBus.Ionomer\">Ionomer</a>, and
  <a href=\"modelica://FCSys.Conditions.ByConnector.ChemicalBus.Liquid\">Liquid</a> models
  are instances of the <a href=\"modelica://FCSys.Conditions.ByConnector.Chemical.Species\">Conditions.ByConnector.Chemical.Species</a>
  model rather than <a href=\"modelica://FCSys.Conditions.ByConnector.Chemical.Reaction\">Conditions.ByConnector.Chemical.Reaction</a>).
  That means that the subconnectors in the
  (<code>chemical</code> connectors of the models in this package are
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors
  (rather than <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a>).</p></html>"),
          Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={2,157,21},
              fillPattern=FillPattern.Solid,
              fillColor={38,196,52},
              lineThickness=0.5)}));

    end PhysicalBus;

    package Physical
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Physical\">Physical</a> connector</html>"
      extends Modelica.Icons.Package;

      model Potential "Specify chemical potential (measure current)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Potential,
          u(final unit="l2.m/(N.T2)"),
          final y(final unit="N/T") = physical.Ndot);

      equation
        physical.mu = u_final;
        annotation (defaultComponentPrefixes="replaceable");
      end Potential;

      model Current "Specify current (measure chemical potential)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Current,
          u(final unit="N/T"),
          final y(final unit="l2.m/(N.T2)") = physical.mu);

      equation
        physical.Ndot = u_final;
        annotation (defaultComponentPrefixes="replaceable");
      end Current;

      model Custom "Custom"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=physical.Ndot);

        Real x=physical.mu "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="physical",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>physical.mu</code> and/or <code>physical.Ndot</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model of a material condition"
          import FCSys.BaseClasses.Utilities.countTrue;
          import FCSys.BaseClasses.Utilities.index;
          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter String formula(start="") "Chemical formula of the species"
            annotation (Dialog(group="Material properties"));
          // The start value prevents a warning in Dymola 7.4.

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification of material condition"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification of material condition",enable=internal),

            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));

          // Properties upon outflow
          parameter Q.Velocity phi[Axis]={0,0,0}
            "<html>Velocity (<b>&phi;</b>)</html>"
            annotation (Dialog(group="Properties upon outflow"));
          parameter Q.PotentialAbsolute sT(start=3000*U.K)
            "Specific entropy-temperature product"
            annotation (Dialog(group="Properties upon outflow"));

          // Included components of translational momentum
          parameter Boolean inclTransX=true "X" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransY=true "Y" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransZ=true "Z" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with
          // the results.

          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));
          Connectors.Physical physical(final n_trans=n_trans,final formula=
                formula) "Connector for phase change"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        protected
          final parameter Integer n_trans=countTrue({inclTransX,inclTransY,
              inclTransZ}) "Number of components of translational momentum";
          final parameter Integer cartTrans[n_trans]=index({inclTransX,
              inclTransY,inclTransZ})
            "Cartesian-axis indices of the components of translational momentum";

          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          physical.phi = phi[cartTrans];
          physical.sT = sT;

          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                  5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (defaultComponentName="physical");
        end PartialCondition;

        type ConditionType = enumeration(
            Current "Specify current (measure chemical potential)",
            Potential "Specify chemical potential (measure current)",
            Custom "Custom") "Types of conditions";

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={2,157,21},
              fillPattern=FillPattern.Solid,
              fillColor={38,196,52})}));

    end Physical;

    package FaceBus
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
      extends Modelica.Icons.Package;

      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"

        extends Modelica.Icons.Package;

        model FaceBus
          "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors, with efforts by default</html>"

          extends FCSys.BaseClasses.Icons.Conditions.Single;

          Phases.Gas gas "Gas" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Graphite graphite "Graphite" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Ionomer ionomer "Ionomer" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Liquid liquid "Liquid" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Connectors.FaceBus negative
            "Negative-side multi-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-100,0})));
          Connectors.FaceBus positive
            "Positive-side multi-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={100,0})));
          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                origin={0,50},
                rotation=270)));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                origin={0,-50},
                rotation=270)));

        equation
          // Gas
          connect(gas.negative, negative.gas) annotation (Line(
              points={{-10,6.10623e-16},{-100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(gas.positive, positive.gas) annotation (Line(
              points={{10,6.10623e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(u.gas, gas.u) annotation (Line(
              points={{5.55112e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(gas.y, y.gas) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));

          // Graphite
          connect(graphite.negative, negative.graphite) annotation (Line(
              points={{-10,6.10623e-16},{-100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(graphite.positive, positive.graphite) annotation (Line(
              points={{10,6.10623e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(u.graphite, graphite.u) annotation (Line(
              points={{5.55112e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(graphite.y, y.graphite) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));

          // Ionomer
          connect(ionomer.negative, negative.ionomer) annotation (Line(
              points={{-10,6.10623e-16},{-100,5.55112e-16}},
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
          connect(u.ionomer, ionomer.u) annotation (Line(
              points={{5.55112e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(ionomer.y, y.ionomer) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));

          // Liquid
          connect(liquid.negative, negative.liquid) annotation (Line(
              points={{-10,6.10623e-16},{-100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(liquid.positive, positive.liquid) annotation (Line(
              points={{10,6.10623e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));
          connect(u.liquid, liquid.u) annotation (Line(
              points={{5.55112e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(liquid.y, y.liquid) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (defaultComponentName="face");
        end FaceBus;

        model FaceBusIsolated
          "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors, with zero normal velocity and otherwise zero flows by default</html>"

          extends FaceBus(
            gas(
              H2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
              H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
              N2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
              O2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal)),
            graphite('C+'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
                'e-'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal)),
            ionomer(
              'C19HF37O5S-'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
              'H+'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal),
              H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal)),
            liquid(H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Material.Current material,
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal)));
          // See note in ChemicalReactionNoFlow.
          annotation (defaultComponentName="face");

        end FaceBusIsolated;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
          extends Modelica.Icons.Package;

          model Gas "Condition for gas"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2=false
              "<html>Hydrogen (H<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face H2 if inclH2 "<html>H<sub>2</sub> conditions</html>"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclN2=false
              "<html>Nitrogen (N<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face N2 if inclN2 "<html>N<sub>2</sub> conditions</html>"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclN2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclO2=false
              "<html>Oxygen (O<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face O2 if inclO2 "<html>O<sub>2</sub> conditions</html>"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclO2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2
            connect(H2.negative, negative.H2) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(H2.positive, positive.H2) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2, H2.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2.y, y.H2) annotation (Line(
                points={{6.10623e-16,-5},{-4.87687e-22,-50},{-4.87687e-22,-50},
                    {5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // H2O
            connect(H2O.negative, negative.H2O) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(H2O.positive, positive.H2O) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2O, H2O.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{6.10623e-16,-5},{-4.87687e-22,-50},{-4.87687e-22,-50},
                    {5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // N2
            connect(N2.negative, negative.N2) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(N2.positive, positive.N2) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.N2, N2.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(N2.y, y.N2) annotation (Line(
                points={{6.10623e-16,-5},{-4.87687e-22,-50},{-4.87687e-22,-50},
                    {5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // O2
            connect(O2.negative, negative.O2) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(O2.positive, positive.O2) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.O2, O2.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(O2.y, y.O2) annotation (Line(
                points={{6.10623e-16,-5},{5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));
            annotation (Diagram(graphics));
          end Gas;

          model Graphite "Condition for graphite"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclC+'=false
              "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face 'C+' if 'inclC+'
              "<html>C<sup>+</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclC+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean 'incle-'=false
              "<html>Electrons (e<sup>-</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face 'e-' if 'incle-'
              "<html>e<sup>-</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='incle-'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

          equation
            // C+
            connect('C+'.negative, negative.'C+') annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('C+'.positive, positive.'C+') annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'C+', 'C+'.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('C+'.y, y.'C+') annotation (Line(
                points={{6.10623e-16,-5},{-4.87687e-22,-50},{-4.87687e-22,-50},
                    {5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // e-
            connect('e-'.negative, negative.'e-') annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('e-'.positive, positive.'e-') annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'e-', 'e-'.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('e-'.y, y.'e-') annotation (Line(
                points={{6.10623e-16,-5},{-4.87687e-22,-50},{-4.87687e-22,-50},
                    {5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Graphite;

          model Ionomer "Condition for ionomer"

            extends BaseClasses.EmptyPhase;

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

            Face.Pair.Face 'C19HF37O5S-' if 'inclC19HF37O5S-'
              "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> conditions</html>"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclC19HF37O5S-'), Placement(transformation(extent={{-10,
                      -10},{10,10}})));

            parameter Boolean 'inclH+'=false
              "<html>Protons (H<sup>+</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face 'H+' if 'inclH+'
              "<html>H<sup>+</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclH+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // C19HF37O5S-
            connect('C19HF37O5S-'.negative, negative.'C19HF37O5S-') annotation
              (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('C19HF37O5S-'.positive, positive.'C19HF37O5S-') annotation
              (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('C19HF37O5S-'.y, y.'C19HF37O5S-') annotation (Line(
                points={{6.10623e-16,-5},{5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // H+
            connect('H+'.negative, negative.'H+') annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('H+'.positive, positive.'H+') annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'H+', 'H+'.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('H+'.y, y.'H+') annotation (Line(
                points={{6.10623e-16,-5},{5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // H2O
            connect(H2O.negative, negative.H2O) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(H2O.positive, positive.H2O) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2O, H2O.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{6.10623e-16,-5},{5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Ionomer;

          model Liquid "Condition for liquid"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2O
            connect(H2O.negative, negative.H2O) annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(H2O.positive, positive.H2O) annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2O, H2O.u) annotation (Line(
                points={{5.55112e-16,50},{5.55112e-16,4},{6.10623e-16,4},{
                    6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{6.10623e-16,-5},{6.10623e-16,-4},{5.55112e-16,-4},{
                    5.55112e-16,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Liquid;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            model EmptyPhase "Empty condition for a phase (no species)"
              extends FCSys.BaseClasses.Icons.Conditions.Pair;

              Connectors.FaceBus negative
                "Negative-side multi-species connector for material, momentum, and energy"
                annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                    iconTransformation(extent={{-110,-10},{-90,10}})));
              Connectors.FaceBus positive
                "Positive-side multi-species connector for material, momentum, and energy"
                annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                    iconTransformation(extent={{90,-10},{110,10}})));
              Connectors.RealInputBus u
                "Input bus for values of specified conditions" annotation (
                  Placement(transformation(
                    extent={{-10,-10},{10,10}},
                    rotation=270,
                    origin={0,50}), iconTransformation(
                    extent={{-10,-10},{10,10}},
                    rotation=270,
                    origin={0,50})));

              Connectors.RealOutputBus y "Output bus of measurements"
                annotation (Placement(transformation(
                    extent={{-10,-10},{10,10}},
                    rotation=270,
                    origin={0,-50}),iconTransformation(
                    extent={{-10,-10},{10,10}},
                    rotation=270,
                    origin={0,-50})));
              annotation (Icon(graphics));

            end EmptyPhase;

          end BaseClasses;

        end Phases;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"

        extends Modelica.Icons.Package;

        model FaceBus
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with efforts by default</html>"

          extends FCSys.BaseClasses.Icons.Conditions.Single;

          Phases.Gas gas "Gas" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Graphite graphite "Graphite" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Ionomer ionomer "Ionomer" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Liquid liquid "Liquid" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Connectors.FaceBus face
            "Connector for material, momentum, and energy of multiple species"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
                iconTransformation(extent={{-10,-50},{10,-30}})));
          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));

        equation
          // Gas
          connect(gas.face, face.gas) annotation (Line(
              points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));

          connect(u.gas, gas.u) annotation (Line(
              points={{-110,5.55112e-16},{-11,5.55112e-16},{-11,6.10623e-16}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(gas.y, y.gas) annotation (Line(
              points={{11,6.10623e-16},{110,6.10623e-16},{110,5.55112e-16}},
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
              points={{-110,5.55112e-16},{-11,5.55112e-16},{-11,6.10623e-16}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(graphite.y, y.graphite) annotation (Line(
              points={{11,6.10623e-16},{110,6.10623e-16},{110,5.55112e-16}},
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
              points={{-110,5.55112e-16},{-11,5.55112e-16},{-11,6.10623e-16}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(ionomer.y, y.ionomer) annotation (Line(
              points={{11,6.10623e-16},{110,6.10623e-16},{110,5.55112e-16}},
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
              points={{-110,5.55112e-16},{-11,5.55112e-16},{-11,6.10623e-16}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(liquid.y, y.liquid) annotation (Line(
              points={{11,6.10623e-16},{110,6.10623e-16},{110,5.55112e-16}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (defaultComponentName="face");
        end FaceBus;

        model FaceBusIsolated
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with zero normal velocity and otherwise zero flows by default</html>"

          extends FaceBus(
            gas(
              H2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),

              H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),

              N2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),

              O2(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal)),

            graphite('C+'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),
                'e-'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal)),

            ionomer(
              'C19HF37O5S-'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),

              'H+'(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal),

              H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal)),

            liquid(H2O(
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Material.Current material,

                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  following(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Translational.Force
                  preceding(final orientation),
                redeclare replaceable
                  Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal)));

          // See note in ChemicalReactionNoFlow.
          annotation (defaultComponentName="face");

        end FaceBusIsolated;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
          extends Modelica.Icons.Package;

          model Gas "Condition for gas"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2=false
              "<html>Hydrogen (H<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face H2 if inclH2
              "<html>H<sub>2</sub> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclN2=false
              "<html>Nitrogen (N<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face N2 if inclN2
              "<html>N<sub>2</sub>Conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclN2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclO2=false
              "<html>Oxygen (O<sub>2</sub>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face O2 if inclO2
              "<html>O<sub>2</sub> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclO2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2
            connect(H2.face, face.H2) annotation (Line(
                points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2, H2.u) annotation (Line(
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2.y, y.H2) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(N2.y, y.N2) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(O2.y, y.O2) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));
            annotation (Diagram(graphics));
          end Gas;

          model Graphite "Condition for graphite"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclC+'=false
              "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face 'C+' if 'inclC+'
              "<html>C<sup>+</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclC+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean 'incle-'=false
              "<html>Electrons (e<sup>-</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face 'e-' if 'incle-'
              "<html>e<sup>-</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='incle-'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

          equation
            // C+
            connect('C+'.face, face.'C+') annotation (Line(
                points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'C+', 'C+'.u) annotation (Line(
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('C+'.y, y.'C+') annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('e-'.y, y.'e-') annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Graphite;

          model Ionomer "Condition for ionomer"

            extends BaseClasses.EmptyPhase;

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

            Face.Single.Face 'C19HF37O5S-' if 'inclC19HF37O5S-'
              "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> conditions</html>"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclC19HF37O5S-'), Placement(transformation(extent={{-10,
                      -10},{10,10}})));

            parameter Boolean 'inclH+'=false
              "<html>Protons (H<sup>+</sup>)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face 'H+' if 'inclH+'
              "<html>H<sup>+</sup> conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable='inclH+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // C19HF37O5S-
            connect('C19HF37O5S-'.face, face.'C19HF37O5S-') annotation (Line(
                points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('C19HF37O5S-'.y, y.'C19HF37O5S-') annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('H+'.y, y.'H+') annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
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
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Ionomer;

          model Liquid "Condition for liquid"

            extends BaseClasses.EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2O=false
              "<html>Water (H<sub>2</sub>O)</html>" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Face H2O if inclH2O
              "<html>H<sub>2</sub>O conditions</html>" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2O
            connect(H2O.face, face.H2O) annotation (Line(
                points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));
            connect(u.H2O, H2O.u) annotation (Line(
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{11,6.10623e-16},{11,0},{100,0},{100,5.55112e-16}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

          end Liquid;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            model EmptyPhase "Empty condition for a phase (no species)"
              extends FCSys.BaseClasses.Icons.Conditions.Single;

              Connectors.FaceBus face
                "Multi-species connector for material, momentum, and energy"
                annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
              Connectors.RealInputBus u
                "Input bus for values of specified conditions" annotation (
                  Placement(transformation(
                    extent={{-10,-10},{10,10}},
                    rotation=0,
                    origin={-100,0}), iconTransformation(
                    extent={{-10,-10},{10,10}},
                    rotation=0,
                    origin={-110,0})));

              Connectors.RealOutputBus y "Output bus of measurements"
                annotation (Placement(transformation(
                    extent={{-10,-10},{10,10}},
                    rotation=0,
                    origin={100,0}),iconTransformation(
                    extent={{-10,-10},{10,10}},
                    rotation=0,
                    origin={110,0})));
              annotation (Icon(graphics));

            end EmptyPhase;

          end BaseClasses;

        end Phases;

      end Single;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={191,191,191},
              lineThickness=0.5)}));

    end FaceBus;

    package Face
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"
      extends Modelica.Icons.Package;
      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors</html>"
        extends Modelica.Icons.Package;
        model Face
          "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors</html>"

          extends FCSys.BaseClasses.Icons.Conditions.Pair;

          replaceable Material.Density material constrainedby
            Conditions.ByConnector.Face.Pair.Material.BaseClasses.PartialCondition
            "Material" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-50,14},{-30,34}})));
          replaceable Translational.Velocity normal(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Pair.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.normal) "Normal translational"
            annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-30,2},{-10,22}})));
          replaceable Translational.Velocity following(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Pair.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.following)
            "<html>1<sup>st</sup> transverse</html>" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-10,-10},{10,10}})));
          replaceable Translational.Velocity preceding(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Pair.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.following)
            "<html>2<sup>nd</sup> transverse</html>" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{10,-22},{30,-2}})));
          replaceable Thermal.Temperature thermal constrainedby
            Conditions.ByConnector.Face.Pair.Thermal.BaseClasses.PartialCondition
            "Thermal" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{30,-34},{50,-14}})));
          // Note:  In Dymola 7.4, the value of y must be specified here instead
          // of at the lower level (e.g., Thermal.Temperature) so that the source
          // subcomponent can be replaced by blocks that don't contain the
          // parameter y.

          Connectors.Face negative
            "Negative-side single-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=180,
                origin={-100,0})));

          Connectors.Face positive
            "Positive-side single-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{90,-10},{110,10}})));

          Connectors.RealInputBus u
            "Input bus for values of specified conditions" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50})));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50})));

        equation
          // Material
          connect(negative, material.negative) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,24},{-50,24}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(material.positive, positive) annotation (Line(
              points={{-30,24},{80,24},{80,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.material, material.u) annotation (Line(
              points={{5.55112e-16,50},{5.55112e-16,30},{-40,30},{-40,29}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(material.y, y.material) annotation (Line(
              points={{-40,19},{-40,-30},{5.55112e-16,-30},{5.55112e-16,-50}},
              color={0,0,127},
              smooth=Smooth.None));

          // Normal
          connect(negative, normal.negative) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,12},{-30,12}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(normal.positive, positive) annotation (Line(
              points={{-10,12},{80,12},{80,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.normal, normal.u) annotation (Line(
              points={{5.55112e-16,50},{5.55112e-16,30},{-20,30},{-20,17}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(normal.y, y.normal) annotation (Line(
              points={{-20,7},{-20,-30},{5.55112e-16,-30},{5.55112e-16,-50}},
              color={0,0,127},
              smooth=Smooth.None));

          // 1st transverse
          connect(negative, following.negative) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,6.10623e-16},{-10,
                  6.10623e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(following.positive, positive) annotation (Line(
              points={{10,6.10623e-16},{80,6.10623e-16},{80,5.55112e-16},{100,
                  5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(u.following, following.u) annotation (Line(
              points={{5.55112e-16,50},{5.55112e-16,30},{6.10623e-16,30},{
                  6.10623e-16,5}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(following.y, y.following) annotation (Line(
              points={{6.10623e-16,-5},{6.10623e-16,-30},{5.55112e-16,-30},{
                  5.55112e-16,-50}},
              color={0,0,127},
              smooth=Smooth.None));

          // 2nd transverse
          connect(negative, preceding.negative) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-12},{10,-12}},

              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(preceding.positive, positive) annotation (Line(
              points={{30,-12},{80,-12},{80,0},{100,0},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(u.preceding, preceding.u) annotation (Line(
              points={{5.55112e-16,50},{5.55112e-16,30},{20,30},{20,-7}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(preceding.y, y.preceding) annotation (Line(
              points={{20,-17},{20,-30},{5.55112e-16,-30},{5.55112e-16,-50}},
              color={0,0,127},
              smooth=Smooth.None));

          // Thermal
          connect(negative, thermal.negative) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-24},{30,-24}},

              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(thermal.positive, positive) annotation (Line(
              points={{50,-24},{80,-24},{80,0},{100,0},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.thermal, thermal.u) annotation (Line(
              points={{5.55112e-16,50},{5.55112e-16,30},{40,30},{40,-19}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(thermal.y, y.thermal) annotation (Line(
              points={{40,-29},{40,-30},{5.55112e-16,-30},{5.55112e-16,-50}},
              color={0,0,127},
              smooth=Smooth.None));

        end Face;

        model FaceIsolated
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector, with zero normal velocity and otherwise zero flows by default</html>"

          extends Face(
            redeclare replaceable
              Conditions.ByConnector.Face.Pair.Material.Current material,
            redeclare replaceable
              Conditions.ByConnector.Face.Pair.Translational.Force following(
                final orientation),
            redeclare replaceable
              Conditions.ByConnector.Face.Pair.Translational.Force preceding(
                final orientation),
            redeclare replaceable
              Conditions.ByConnector.Face.Pair.Thermal.HeatRate thermal);
          // See note in ChemicalReactionNoFlow.
          annotation (defaultComponentName="face");

        end FaceIsolated;

        package Material "Material conditions"
          extends Modelica.Icons.Package;

          model Density
            "Specify density difference (measure diffusion current), with conservation of material"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Density,
              u(final unit="N/l3"),
              final y(final unit="N/T") = negative.Ndot);

          equation
            positive.rho - negative.rho = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="material");
          end Density;

          model Pressure
            "Specify pressure difference (measure diffusion current)"
            extends Single.Material.BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Pressure,
              u(final unit="m/(l.T2)"),
              final y(final unit="N/T") = face.Ndot);

            replaceable package Data =
                Characteristics.BaseClasses.Characteristic constrainedby
              Characteristics.BaseClasses.CharacteristicEOS
              "Characteristic data" annotation (
              Dialog(group="Material properties"),
              choicesAllMatching=true,
              __Dymola_choicesFromPackage=true,
              Placement(transformation(extent={{-60,40},{-40,60}}),
                  iconTransformation(extent={{-10,90},{10,110}})));

          equation
            Data.p_Tv(positive.T, 1/positive.rho) - Data.p_Tv(negative.T, 1/
              negative.rho) = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="material");
          end Pressure;

          model Current
            "Specify diffusion current (measure density difference), with conservation of material"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Current,
              u(final unit="N/T"),
              final y(final unit="N/l3") = positive.rho - negative.rho);

          equation
            negative.Ndot = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="material");
          end Current;

          model Custom "Custom"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=negative.Ndot);

            Real x=positive.rho - negative.rho
              "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="material",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>negative.rho</code>, <code>positive.rho</code>, <code>negative.Ndot</code>,
    and/or <code>positive.Ndot</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a normal condition"

              extends Pair.BaseClasses.PartialCondition;

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // Conservation of material
              0 = negative.Ndot + positive.Ndot;

              // No flows of other quantities
              // ----------------------------
              // Translational momentum
              negative.mPhidot = {0,0,0};
              positive.mPhidot = {0,0,0};
              //
              // Thermal energy
              negative.Qdot = 0;
              positive.Qdot = 0;
              annotation (defaultComponentName="material");
            end PartialCondition;

            type ConditionType = enumeration(
                Density "Specify density difference (measure current)",
                Pressure "Specify pressure difference (measure current)",
                Current "Specify current (measure density difference)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Material;

        package Translational "Translational conditions"
          extends Modelica.Icons.Package;

          model Velocity
            "Specify velocity difference (measure force), with conservation of translational momentum"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Velocity,
              u(final unit="l/T"),
              final y(final unit="l.m/T2") = negative.mPhidot[orientation]);

          equation
            positive.phi[orientation] - negative.phi[orientation] = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="translational");
          end Velocity;

          model Force
            "Specify force (measure velocity difference), with conservation of translational momentum"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Force,
              u(final unit="l.m/T2"),
              final y(final unit="l/T") = positive.phi[orientation] - negative.phi[
                orientation]);

          equation
            negative.mPhidot[orientation] = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="translational");
          end Force;

          model Custom "Custom"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=negative.mPhidot[
                  orientation]);

            Real x=positive.phi[orientation] - negative.phi[orientation]
              "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="translational",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>negative.phi[orientation]</code>, <code>positive.phi[orientation]</code>,
    <code>negative.mPhidot[orientation]</code> and/or <code>positive.mPhidot[orientation]</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a transverse translational condition"
              import FCSys.BaseClasses.Utilities.cartWrap;

              extends Pair.BaseClasses.PartialCondition;

              parameter Orientation orientation=Orientation.normal
                "Orientation of translational momentum";

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // Conservation of translational momentum in the present direction
              0 = negative.mPhidot[orientation] + positive.mPhidot[orientation];

              // No flows of other quantities
              // ----------------------------
              // Material
              negative.Ndot = 0;
              positive.Ndot = 0;
              //
              // Translational momentum in the other directions
              negative.mPhidot[cartWrap(orientation + 1)] = 0;
              positive.mPhidot[cartWrap(orientation + 1)] = 0;
              negative.mPhidot[cartWrap(orientation - 1)] = 0;
              positive.mPhidot[cartWrap(orientation - 1)] = 0;
              //
              // Heat
              negative.Qdot = 0;
              positive.Qdot = 0;
              annotation (defaultComponentName="translational");
            end PartialCondition;

            type ConditionType = enumeration(
                Velocity "Specify velocity difference (measure force)",
                Force "Specify force (measure velocity difference)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Translational;

        package Thermal "Thermal conditions"
          extends Modelica.Icons.Package;

          model Temperature
            "Specify temperature difference (measure heat flow rate), with conservation of energy"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Temperature,
              u(final unit="l2.m/(N.T2)", displayUnit="K"),
              final y(unit="l2.m/T3") = negative.Qdot);

          equation
            positive.T - negative.T = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="thermal");
          end Temperature;

          model HeatRate
            "Specify heat flow rate (measure temperature difference), with conservation of energy"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.HeatRate,
              u(final unit="l2.m/T3"),
              final y(
                final unit="l2.m/(N.T2)",
                displayUnit="K") = positive.T - negative.T);

          equation
            negative.Qdot = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="thermal");
          end HeatRate;

          model Custom
            "Apply condition to a custom expression, with conservation of energy"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=negative.Qdot);

            Real x=positive.T - negative.T
              "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>negative.T</code>, <code>positive.T</code>, <code>negative.Qdot</code> and/or <code>positive.Qdot</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a thermal condition"

              extends Pair.BaseClasses.PartialCondition;

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // Conservation of energy (no storage)
              0 = negative.Qdot + positive.Qdot;

              // No flows of other quantities
              // ----------------------------
              // Material
              negative.Ndot = 0;
              positive.Ndot = 0;
              //
              // Translational momentum
              negative.mPhidot = {0,0,0};
              positive.mPhidot = {0,0,0};
              annotation (defaultComponentName="thermal");
            end PartialCondition;

            type ConditionType = enumeration(
                Temperature
                  "Specify temperature difference (measure heat flow rate)",
                HeatRate
                  "Specify heat flow rate (measure temperature difference)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Thermal;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;

          partial model PartialCondition
            "Partial model to specify and measure conditions on a pair of connectors"
            extends FCSys.BaseClasses.Icons.Conditions.Pair;

            parameter Boolean internal=true "Use internal specification"
              annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(group="Specification"));

            replaceable Modelica.Blocks.Sources.RealExpression source if
              internal constrainedby Modelica.Blocks.Interfaces.SO
              "Source of internal specification" annotation (
              __Dymola_choicesFromPackage=true,
              Dialog(group="Specification",enable=internal),
              Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={40,20})));

            Connectors.RealInput u if not internal
              "Value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,50})));

            Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                  group="Measurement"), Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,-50}), iconTransformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,-50})));

            Connectors.Face negative
              "Negative-side connector to transport material, momentum, and energy of a single species"
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
            Connectors.Face positive
              "Positive-side connector to transport material, momentum, and energy of a single species"
              annotation (Placement(transformation(extent={{90,-10},{110,10}})));

          protected
            Connectors.RealOutputInternal u_final
              "Final value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,-20})));

          equation
            connect(u, u_final) annotation (Line(
                points={{5.55112e-16,50},{0,0},{0,-20},{5.55112e-16,-20}},
                color={0,0,127},
                smooth=Smooth.None));

            connect(source.y, u_final) annotation (Line(
                points={{40,9},{40,0},{0,0},{0,-20},{5.55112e-16,-20}},
                color={0,0,127},
                smooth=Smooth.None));
            annotation (Icon(graphics));
          end PartialCondition;

        end BaseClasses;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"
        extends Modelica.Icons.Package;
        model Face
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"

          extends FCSys.BaseClasses.Icons.Conditions.Single;

          replaceable Material.Density material(source(y=4*U.C/U.cc))
            constrainedby
            Conditions.ByConnector.Face.Single.Material.BaseClasses.PartialCondition
            "Material" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-66,22},{-46,42}})));
          replaceable Translational.Velocity normal(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Single.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.normal) "Normal translational"
            annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-38,10},{-18,30}})));
          replaceable Translational.Velocity following(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Single.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.following)
            "<html>1<sup>st</sup> transverse</html>" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{-10,-2},{10,18}})));
          replaceable Translational.Velocity preceding(final orientation)
            constrainedby
            Conditions.ByConnector.Face.Single.Translational.BaseClasses.PartialCondition(
              orientation=Orientation.preceding)
            "<html>2<sup>nd</sup> transverse</html>" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{18,-18},{38,2}})));
          replaceable Thermal.Temperature thermal(source(y=298.15*U.K))
            constrainedby
            Conditions.ByConnector.Face.Single.Thermal.BaseClasses.PartialCondition
            "Thermal" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Conditions"),
            Placement(transformation(extent={{46,-30},{66,-10}})));
          // Note:  In Dymola 7.4, the value of y must be specified here instead
          // of at the lower level (e.g., Thermal.Temperature) so that the source
          // subcomponent can be replaced by blocks that don't contain the
          // parameter y.

          Connectors.Face face
            "Single-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
                iconTransformation(extent={{-10,-50},{10,-30}})));

          Connectors.RealInputBus u
            "Input bus for values of specified conditions" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-100,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={100,0}),iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));

        equation
          // Material
          connect(material.face, face) annotation (Line(
              points={{-56,28},{-56,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.material, material.u) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,32},{-67,32}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%first",
              index=-1,
              extent={{-1,3},{-1,3}}));
          connect(material.y, y.material) annotation (Line(
              points={{-45,32},{80,32},{80,0},{100,0},{100,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{1,3},{1,3}}));

          // Normal translational
          connect(normal.face, face) annotation (Line(
              points={{-28,16},{-28,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.normal, normal.u) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,20},{-39,20}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%first",
              index=-1,
              extent={{-1,3},{-1,3}}));
          connect(normal.y, y.normal) annotation (Line(
              points={{-17,20},{80,20},{80,5.55112e-16},{100,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{1,3},{1,3}}));

          // 1st transverse
          connect(following.face, face) annotation (Line(
              points={{6.10623e-16,4},{6.10623e-16,-30},{0,-30},{0,-40},{
                  5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(u.following, following.u) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,8},{-11,8}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%first",
              index=-1,
              extent={{-1,3},{-1,3}}));
          connect(following.y, y.following) annotation (Line(
              points={{11,8},{80,8},{80,5.55112e-16},{100,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{1,3},{1,3}}));

          // 2nd transverse
          connect(preceding.face, face) annotation (Line(
              points={{28,-12},{28,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.preceding, preceding.u) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-8},{17,-8}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%first",
              index=-1,
              extent={{-1,3},{-1,3}}));
          connect(preceding.y, y.preceding) annotation (Line(
              points={{39,-8},{80,-8},{80,0},{100,0},{100,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{1,3},{1,3}}));

          // Thermal
          connect(thermal.face, face) annotation (Line(
              points={{56,-24},{56,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
              color={127,127,127},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(u.thermal, thermal.u) annotation (Line(
              points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-20},{45,-20}},

              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%first",
              index=-1,
              extent={{-1,3},{-1,3}}));

          connect(thermal.y, y.thermal) annotation (Line(
              points={{67,-20},{80,-20},{80,0},{100,0},{100,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{1,3},{1,3}}));

        end Face;

        model FaceIsolated
          extends Face(
            redeclare replaceable
              Conditions.ByConnector.Face.Single.Material.Current material,
            redeclare replaceable
              Conditions.ByConnector.Face.Single.Translational.Force following(
                final orientation),
            redeclare replaceable
              Conditions.ByConnector.Face.Single.Translational.Force preceding(
                final orientation),
            redeclare replaceable
              Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal);
          // See note in ChemicalReactionNoFlow.
          annotation (defaultComponentName="face");

        end FaceIsolated;

        package Material "Material conditions"
          extends Modelica.Icons.Package;

          model Density "Specify density (measure current)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Density,
              u(final unit="N/l3"),
              final y(final unit="N/T") = face.Ndot);

          equation
            face.rho = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="material");
          end Density;

          model Pressure "Specify pressure (measure diffusion current)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Pressure,
              u(final unit="m/(l.T2)"),
              final y(final unit="N/T") = face.Ndot);

            replaceable package Data =
                Characteristics.BaseClasses.Characteristic constrainedby
              Characteristics.BaseClasses.Characteristic "Characteristic data"
              annotation (
              Dialog(group="Material properties"),
              choicesAllMatching=true,
              __Dymola_choicesFromPackage=true,
              Placement(transformation(extent={{-60,40},{-40,60}}),
                  iconTransformation(extent={{-10,90},{10,110}})));

          equation
            Data.p_Tv(face.T, 1/face.rho) = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="material",
              Documentation(info="<html>
  <p>The default characteristic data represents an ideal gas.</p></html>"));
          end Pressure;

          model Current "Specify diffusion current (measure density)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Current,
              u(final unit="N/T"),
              final y(final unit="N/l3") = face.rho);

          equation
            face.Ndot = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="material");
          end Current;

          model Custom "Custom"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=face.Ndot);

            Real x=face.rho "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="material",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.rho</code> and/or <code>face.Ndot</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a normal condition"

              extends Single.BaseClasses.PartialCondition;

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // No flows of other quantities
              face.mPhidot = {0,0,0} "Translational momentum";
              face.Qdot = 0 "Heat";
              annotation (defaultComponentName="material");
            end PartialCondition;

            type ConditionType = enumeration(
                Density "Specify density (measure current)",
                Pressure "Specify pressure (measure current)",
                Current "Specify current (measure Denis)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Material;

        package Translational "Translational conditions"
          extends Modelica.Icons.Package;

          model Velocity "Specify velocity (measure force)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Velocity,
              u(final unit="l/T"),
              final y(final unit="l.m/T2") = face.mPhidot[orientation]);

          equation
            face.phi[orientation] = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="translational");
          end Velocity;

          model Current "Specify advective current (measure force)"
            import assert = FCSys.BaseClasses.Utilities.assertEval;
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Current,
              u(final unit="N/T"),
              final y(final unit="l.m/T2") = face.mPhidot[orientation]);

            Q.Area A=U.cm^2 "Cross-sectional area";

          initial equation
            assert(orientation == Orientation.normal,
              "FCSys.Conditions.ByConnector.Face.Single.Translational.Current is only intended for the normal direction.");

          equation
            face.phi[orientation]*face.rho*A = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="translational",
              Documentation(info="<html><p>The advective current is in the globally 
    positive direction (not into the component).
    This model is meaningful only for the normal direction
  (<code>orientation = Orientation.normal).</p></html>"));
          end Current;

          model Force "Specify force (measure velocity)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Force,
              u(final unit="l.m/T2"),
              final y(final unit="l/T") = face.phi[orientation]);

          equation
            face.mPhidot[orientation] = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="translational");
          end Force;

          model Custom "Custom"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=face.mPhidot[orientation]);

            Real x=face.phi[orientation]
              "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="translational",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.phi[orientation]</code> and/or <code>face.mPhidot[orientation]</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a transverse translational condition"
              import FCSys.BaseClasses.Utilities.cartWrap;

              extends Single.BaseClasses.PartialCondition;

              parameter Orientation orientation=Orientation.normal
                "Orientation of translational momentum";

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // No flows of other quantities
              face.Ndot = 0 "Material";
              face.mPhidot[cartWrap(orientation + 1)] = 0
                "Translational momentum in the following direction";
              face.mPhidot[cartWrap(orientation - 1)] = 0
                "Translational momentum in the preceding direction";
              face.Qdot = 0 "Heat";
              annotation (defaultComponentName="translational");
            end PartialCondition;

            type ConditionType = enumeration(
                Velocity "Specify velocity (measure force)",
                Current "Specify current (measure force)",
                Force "Specify force (measure velocity)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Translational;

        package Thermal "Thermal conditions"
          extends Modelica.Icons.Package;

          model Temperature "Specify temperature (measure heat flow rate)"
            extends BaseClasses.PartialCondition(
              final conditionType=BaseClasses.ConditionType.Temperature,
              u(final unit="l2.m/(N.T2)", displayUnit="K"),
              final y(unit="l2.m/T3") = face.Qdot);

          equation
            face.T = u_final;
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="thermal");
          end Temperature;

          model HeatRate "Specify heat flow rate (measure temperature)"
            extends BaseClasses.PartialCondition(
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

          model Custom "Custom"
            extends BaseClasses.PartialCondition(final conditionType=
                  BaseClasses.ConditionType.Custom, y=face.Qdot);

            Real x=face.T "Expression to which the condition is applied"
              annotation (Dialog(group="Specification"));

          equation
            x = u_final;
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal",
              Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.T</code> and/or <code>face.Qdot</code>.</p></html>"));
          end Custom;

          package BaseClasses "Base classes (generally not for direct use)"
            extends Modelica.Icons.BasesPackage;
            partial model PartialCondition
              "Partial model for a thermal condition"

              extends Single.BaseClasses.PartialCondition;

              constant ConditionType conditionType "Type of condition";
              // Note:  This is included so that the type of condition is recorded with
              // the results.

            equation
              // No flows of other quantities
              face.Ndot = 0 "Material";
              face.mPhidot = {0,0,0} "Translational momentum";
              annotation (defaultComponentName="thermal");
            end PartialCondition;

            type ConditionType = enumeration(
                Temperature "Specify temperature (measure heat flow rate)",
                HeatRate "Specify heat flow rate (measure temperature)",
                Custom "Custom") "Types of conditions";

          end BaseClasses;

        end Thermal;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;

          partial model PartialCondition
            "Partial model to specify and measure conditions on a connector"
            extends FCSys.BaseClasses.Icons.Conditions.Single;

            parameter Boolean internal=true "Use internal specification"
              annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(group="Specification"));

            replaceable Modelica.Blocks.Sources.RealExpression source if
              internal constrainedby Modelica.Blocks.Interfaces.SO
              "Source of internal specification" annotation (
              __Dymola_choicesFromPackage=true,
              Dialog(group="Specification",enable=internal),
              Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-70,30})));

            Connectors.RealInput u if not internal
              "Value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-110,0})));

            Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                  group="Measurement"), Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={110,0}), iconTransformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={110,0})));

            Connectors.Face face
              "Connector to transport material, momentum, and energy of a single species"
              annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

          protected
            Connectors.RealOutputInternal u_final
              "Final value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-20,0})));

          equation
            connect(u, u_final) annotation (Line(
                points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},
                    {-20,5.55112e-16}},
                color={0,0,127},
                smooth=Smooth.None));

            connect(source.y, u_final) annotation (Line(
                points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},

                color={0,0,127},
                smooth=Smooth.None));

            annotation (Icon(graphics));
          end PartialCondition;

        end BaseClasses;

      end Single;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={191,191,191})}));

    end Face;

    package Inert
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector</html>"
      extends Modelica.Icons.Package;

      model Inert
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector, with efforts by default</html>"
        import FCSys.BaseClasses.Utilities.countTrue;
        import FCSys.BaseClasses.Utilities.enumerate;
        import FCSys.BaseClasses.Utilities.index;
        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        // Conditions
        replaceable Translational.Velocity translationalX(final axis) if
          inclTransX constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.x) "X-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransX),
          Placement(transformation(extent={{-58,8},{-38,28}})));

        replaceable Translational.Velocity translationalY(final axis) if
          inclTransY constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.y) "Y-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransY),
          Placement(transformation(extent={{-26,-4},{-6,16}})));

        replaceable Translational.Velocity translationalZ(final axis) if
          inclTransZ constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.z) "Z-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransZ),
          Placement(transformation(extent={{6,-16},{26,4}})));

        replaceable Thermal.Temperature thermal(source(y=298.15*U.K))
          constrainedby
          Conditions.ByConnector.Inert.Thermal.BaseClasses.PartialCondition
          "Thermal" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{38,-30},{58,-10}})));

        Connectors.RealInputBus u
          "Input bus for values of specified conditions" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealOutputBus y "Output bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.Inert inert(final n_trans=countTrue({inclTransX,inclTransY,
              inclTransZ}))
          "Single-species connector for material, momentum, and energy"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-50},{10,-30}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer cartTrans[n_trans]=index({inclTransX,
            inclTransY,inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

      equation
        // X-axis translational
        connect(translationalX.translational, inert.translational) annotation (
            Line(
            points={{-48,14},{-48,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.translationalX, translationalX.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,18},{-59,18}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalX.y, y.translationalX) annotation (Line(
            points={{-37,18},{90,18},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Y-axis translational
        connect(translationalY.translational, inert.translational) annotation (
            Line(
            points={{-16,2},{-16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(u.translationalY, translationalY.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,6},{-27,6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalY.y, y.translationalY) annotation (Line(
            points={{-5,6},{90,6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Z-axis translational
        connect(translationalZ.translational, inert.translational) annotation (
            Line(
            points={{16,-10},{16,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.translationalZ, translationalZ.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,-6},{5,-6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalZ.y, y.translationalZ) annotation (Line(
            points={{27,-6},{90,-6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Thermal
        connect(thermal.thermal, inert.thermal) annotation (Line(
            points={{48,-24},{48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.thermal, thermal.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,-20},{37,-20}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(thermal.y, y.thermal) annotation (Line(
            points={{59,-20},{90,-20},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

      end Inert;

      model InertIsolated
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector, with zero flows by default</html>"

        extends Inert(
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalX(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalY(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalZ(
              final axis),
          redeclare replaceable Conditions.ByConnector.Inert.Thermal.HeatRate
            thermal);
        // See note in ChemicalReactionNoFlow.
        annotation (defaultComponentName="inert");

      end InertIsolated;

      model InertInternal
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector, with efforts by default</html>"
        import FCSys.BaseClasses.Utilities.countTrue;
        import FCSys.BaseClasses.Utilities.enumerate;
        import FCSys.BaseClasses.Utilities.index;
        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        // Included subconnectors
        parameter Boolean inclTranslational=true "Translational" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Included subconnectors",compact=true));
        parameter Boolean inclThermal=true "Thermal" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Included subconnectors",compact=true));

        // Conditions
        replaceable Translational.Velocity translationalX(final axis) if
          inclTransX constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.x) "X-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransX),
          Placement(transformation(extent={{-58,8},{-38,28}})));

        replaceable Translational.Velocity translationalY(final axis) if
          inclTransY constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.y) "Y-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransY),
          Placement(transformation(extent={{-26,-4},{-6,16}})));

        replaceable Translational.Velocity translationalZ(final axis) if
          inclTransZ constrainedby
          Conditions.ByConnector.Inert.Translational.BaseClasses.PartialCondition(
            axis=Axis.z) "Z-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransZ),
          Placement(transformation(extent={{6,-16},{26,4}})));

        replaceable Thermal.Temperature thermal(source(y=298.15*U.K))
          constrainedby
          Conditions.ByConnector.Inert.Thermal.BaseClasses.PartialCondition
          "Thermal" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{38,-30},{58,-10}})));

        Connectors.RealInputBus u
          "Input bus for values of specified conditions" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealOutputBus y "Output bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.InertInternal inert(
          final n_trans=countTrue({inclTransX,inclTransY,inclTransZ}),
          inclTranslational=inclTranslational,
          inclThermal=inclThermal)
          "Single-species connector for material, momentum, and energy"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-50},{10,-30}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer cartTrans[n_trans]=index({inclTransX,
            inclTransY,inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

      equation
        // X-axis translational
        connect(translationalX.translational, inert.translational) annotation (
            Line(
            points={{-48,14},{-48,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.translationalX, translationalX.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,18},{-59,18}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalX.y, y.translationalX) annotation (Line(
            points={{-37,18},{90,18},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Y-axis translational
        connect(translationalY.translational, inert.translational) annotation (
            Line(
            points={{-16,2},{-16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(u.translationalY, translationalY.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,6},{-27,6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalY.y, y.translationalY) annotation (Line(
            points={{-5,6},{90,6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Z-axis translational
        connect(translationalZ.translational, inert.translational) annotation (
            Line(
            points={{16,-10},{16,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.translationalZ, translationalZ.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,-6},{5,-6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalZ.y, y.translationalZ) annotation (Line(
            points={{27,-6},{90,-6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Thermal
        connect(thermal.thermal, inert.thermal) annotation (Line(
            points={{48,-24},{48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(u.thermal, thermal.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-90,0},{-90,-20},{37,-20}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(thermal.y, y.thermal) annotation (Line(
            points={{59,-20},{90,-20},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));
        annotation (defaultComponentName="inert");
      end InertInternal;

      model InertInternalIsolated
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector, with zero flows by default</html>"

        extends InertInternal(
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalX(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalY(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Inert.Translational.Force translationalZ(
              final axis),
          redeclare replaceable Conditions.ByConnector.Inert.Thermal.HeatRate
            thermal);
        // See note in ChemicalReactionNoFlow.
        annotation (defaultComponentName="inert");

      end InertInternalIsolated;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;
        model Velocity "Specify velocity (measure force)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"),
            final y(final unit="l.m/T2") = translational.mPhidot[transCart[axis]]);

        equation
          if n_trans > 0 then
            translational.phi[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Velocity;

        model Force "Specify force (measure velocity)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"),
            final y(final unit="l/T") = translational.phi[transCart[axis]]);

        equation
          if n_trans > 0 then
            translational.mPhidot[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Force;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=translational.mPhidot[transCart[axis]]);

          Real x=translational.phi[transCart[axis]]
            "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          if n_trans > 0 then
            x = u_final;
          end if;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="translational",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>translational.phi[transCart[axis]]</code> and/or <code>translational.mPhidot[transCart[axis]]</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a translational condition"
            extends ByConnector.Inert.BaseClasses.PartialCondition;

            parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

            Connectors.Translational translational(final n_trans=n_trans)
              "Connector to exchange translational momentum" annotation (
                Placement(transformation(extent={{-10,-50},{10,-30}})));

          protected
            outer parameter Integer n_trans
              "Number of components of translational momentum";
            outer parameter Integer cartTrans[:]
              "Cartesian-axis indices of the components of translational momentum";
            outer parameter Integer transCart[Axis]
              "Translational-momentum-component indices of the Cartesian axes";

          equation
            // Zero values of other flows
            for i in 1:n_trans loop
              if cartTrans[i] <> axis then
                translational.mPhidot[i] = 0 "Force along the other axes";
              end if;
            end for;
            annotation (defaultComponentName="translational");
          end PartialCondition;

          type ConditionType = enumeration(
              Velocity "Specify velocity (measure force)",
              Force "Specify force (measure velocity)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Translational;

      package Thermal "Thermal conditions"
        extends Modelica.Icons.Package;

        model Temperature "Specify temperature (measure heat flow rate)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"),
            final y(final unit="l2.m/T3") = thermal.Qdot);

        equation
          thermal.T = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal");
        end Temperature;

        model HeatRate "Specify heat flow rate (measure temperature)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.HeatRate,
            u(final unit="l2.m/T3"),
            final y(
              final unit="l2.m/(N.T2)",
              displayUnit="K") = thermal.T);

        equation
          thermal.Qdot = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal");
        end HeatRate;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=thermal.Qdot);

          Real x=thermal.T "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>thermal.T</code> and/or <code>thermal.Qdot</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a thermal condition"
            extends ByConnector.Inert.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

            Connectors.ThermalDiffusion thermal "Connector to exchange heat"
              annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
            annotation (defaultComponentName="thermal");

          end PartialCondition;

          type ConditionType = enumeration(
              Temperature "Specify temperature (measure heat flow rate)",
              HeatRate "Specify heat flow rate (measure temperature)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Thermal;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model of a condition"

          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification",enable=internal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));
          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));

        protected
          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                  5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));

        end PartialCondition;

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}), Ellipse(
              extent={{-40,20},{20,-40}},
              fillColor={47,107,251},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0})}));

    end Inert;

    package InertDalton
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
      extends Modelica.Icons.Package;

      model InertDalton
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, with efforts by default</html>"
        import FCSys.BaseClasses.Utilities.countTrue;
        import FCSys.BaseClasses.Utilities.enumerate;
        import FCSys.BaseClasses.Utilities.index;
        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        // Conditions
        replaceable Dalton.Volume dalton(source(y=U.cc)) constrainedby
          Conditions.ByConnector.InertDalton.Dalton.BaseClasses.PartialCondition
          "Dalton" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{-74,20},{-54,40}})));
        replaceable Translational.Velocity translationalX(final axis) if
          inclTransX constrainedby
          Conditions.ByConnector.InertDalton.Translational.BaseClasses.PartialCondition(
            axis=Axis.x) "X-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransX),
          Placement(transformation(extent={{-42,8},{-22,28}})));

        replaceable Translational.Velocity translationalY(final axis) if
          inclTransY constrainedby Translational.BaseClasses.PartialCondition(
            axis=Axis.y) "Y-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransY),
          Placement(transformation(extent={{-10,-4},{10,16}})));

        replaceable Translational.Velocity translationalZ(final axis) if
          inclTransZ constrainedby
          Conditions.ByConnector.InertDalton.Translational.BaseClasses.PartialCondition(
            axis=Axis.z) "Z-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransZ),
          Placement(transformation(extent={{22,-16},{42,4}})));

        replaceable Thermal.Temperature thermal(source(y=298.15*U.K))
          constrainedby
          Conditions.ByConnector.InertDalton.Thermal.BaseClasses.PartialCondition
          "Thermal" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions"),
          Placement(transformation(extent={{54,-30},{74,-10}})));

        Connectors.RealInputBus u
          "Input bus for values of specified conditions" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealOutputBus y "Output bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.InertDalton inert(final n_trans=countTrue({inclTransX,
              inclTransY,inclTransZ}))
          "Single-species connector for material, momentum, and energy, with additivity of pressure"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-50},{10,-30}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer cartTrans[n_trans]=index({inclTransX,
            inclTransY,inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

      equation
        // Dalton
        connect(dalton.inert, inert) annotation (Line(
            points={{-64,26},{-64,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={11,43,197},
            smooth=Smooth.None));
        connect(u.dalton, dalton.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,30},{-75,30}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(dalton.y, y.dalton) annotation (Line(
            points={{-53,30},{90,30},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // X-axis translational
        connect(translationalX.inert, inert) annotation (Line(
            points={{-32,14},{-32,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={11,43,197},
            smooth=Smooth.None));
        connect(u.translationalX, translationalX.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,18},{-43,18}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalX.y, y.translationalX) annotation (Line(
            points={{-21,18},{90,18},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Y-axis translational
        connect(translationalY.inert, inert) annotation (Line(
            points={{6.10623e-16,2},{6.10623e-16,-30},{5.55112e-16,-30},{
                5.55112e-16,-40}},
            color={11,43,197},
            smooth=Smooth.None));

        connect(u.translationalY, translationalY.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,6},{-11,6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalY.y, y.translationalY) annotation (Line(
            points={{11,6},{90,6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Z-axis translational
        connect(translationalZ.inert, inert) annotation (Line(
            points={{32,-10},{32,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={11,43,197},
            smooth=Smooth.None));
        connect(u.translationalZ, translationalZ.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,-6},{21,-6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalZ.y, y.translationalZ) annotation (Line(
            points={{43,-6},{90,-6},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Thermal
        connect(thermal.inert, inert) annotation (Line(
            points={{64,-24},{64,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={11,43,197},
            smooth=Smooth.None));
        connect(u.thermal, thermal.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,-20},{53,-20}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(thermal.y, y.thermal) annotation (Line(
            points={{75,-20},{90,-20},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

      end InertDalton;

      model InertDaltonEmpty
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, with zero flows by default</html>"

        extends InertDalton(
          redeclare replaceable
            Conditions.ByConnector.InertDalton.Dalton.Pressure dalton,
          redeclare replaceable
            Conditions.ByConnector.InertDalton.Translational.Force
            translationalX(final axis),
          redeclare replaceable
            Conditions.ByConnector.InertDalton.Translational.Force
            translationalY(final axis),
          redeclare replaceable
            Conditions.ByConnector.InertDalton.Translational.Force
            translationalZ(final axis),
          redeclare replaceable
            Conditions.ByConnector.InertDalton.Thermal.HeatRate thermal);
        // See note in ChemicalReactionNoFlow.
        annotation (defaultComponentName="inertDalton");

      end InertDaltonEmpty;

      package Dalton "Conditions for additivity of volume"
        extends Modelica.Icons.Package;

        model Volume "Specify volume (measure pressure)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Volume,
            u(final unit="l3"),
            final y(final unit="m/(l.T2)") = inert.p);

        equation
          inert.V = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="dalton");
        end Volume;

        model Pressure "Specify pressure (measure volume)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Pressure,
            u(final unit="m/(l.T2)"),
            final y(final unit="l3") = inert.V);

        equation
          inert.p = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="dalton");
        end Pressure;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=inert.p);

          Real x=inert.V "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="dalton",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>inert.V</code> and/or <code>inert.p</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model of a volume/pressure condition"
            extends ByConnector.InertDalton.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          equation
            // Zero values of other flows
            inert.mPhidot = zeros(n_trans) "Force";
            inert.Qdot = 0 "Heat flow rate";
            annotation (defaultComponentName="dalton");
          end PartialCondition;

          type ConditionType = enumeration(
              Volume "Specify volume (measure pressure)",
              Pressure "Specify pressure (measure volume)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Dalton;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;
        model Velocity "Specify velocity (measure force)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"),
            final y(final unit="l.m/T2") = inert.mPhidot[transCart[axis]]);

        equation
          if n_trans > 0 then
            inert.phi[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Velocity;

        model Force "Specify force (measure velocity)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"),
            final y(final unit="l/T") = inert.phi[transCart[axis]]);

        equation
          if n_trans > 0 then
            inert.mPhidot[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Force;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=inert.mPhidot[transCart[axis]]);

          Real x=inert.phi[transCart[axis]]
            "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          if n_trans > 0 then
            x = u_final;
          end if;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="translational",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>inert.phi[transCart[axis]</code> and/or <code>inert.mPhidot[transCart[axis]</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a translational condition"
            extends ByConnector.InertDalton.BaseClasses.PartialCondition;

            parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          protected
            outer parameter Integer cartTrans[:]
              "Cartesian-axis indices of the components of translational momentum";
            outer parameter Integer transCart[Axis]
              "Translational-momentum-component indices of the Cartesian axes";

          equation
            // Zero values of other flows
            inert.p = 0 "Pressure";
            for i in 1:n_trans loop
              if cartTrans[i] <> axis then
                inert.mPhidot[i] = 0 "Force along the other axes";
              end if;
            end for;
            inert.Qdot = 0 "Heat flow rate";
            annotation (defaultComponentName="translational");
          end PartialCondition;

          type ConditionType = enumeration(
              Velocity "Specify velocity (measure force)",
              Force "Specify force (measure velocity)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Translational;

      package Thermal "Thermal conditions"
        extends Modelica.Icons.Package;

        model Temperature "Specify temperature (measure heat flow rate)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"),
            final y(final unit="l2.m/T3") = inert.Qdot);

        equation
          inert.T = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal");
        end Temperature;

        model HeatRate "Specify heat flow rate (measure temperature)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.HeatRate,
            u(final unit="l2.m/T3"),
            final y(
              final unit="l2.m/(N.T2)",
              displayUnit="K") = inert.T);

        equation
          inert.Qdot = u_final;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="thermal");
        end HeatRate;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=inert.Qdot);

          Real x=inert.T "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          x = u_final;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>inert.T</code> and/or <code>inert.Qdot</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a thermal condition"
            extends ByConnector.InertDalton.BaseClasses.PartialCondition;

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          equation
            // Zero values of other flows
            inert.p = 0 "Pressure";
            inert.mPhidot = zeros(n_trans) "Force";
            annotation (defaultComponentName="thermal");
          end PartialCondition;

          type ConditionType = enumeration(
              Temperature "Specify temperature (measure heat flow rate)",
              HeatRate "Specify heat flow rate (measure temperature)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Thermal;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model of a condition"
          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter Boolean inclTransX=true "X" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransY=true "Y" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean inclTransZ=true "Z" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Axes with translational momentum included",
              compact=true));

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification",enable=internal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));
          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));
          Connectors.InertDalton inert(final n_trans=n_trans)
            "Connector for translational momentum and energy, with additivity of pressure"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        protected
          outer parameter Integer n_trans
            "Number of components of translational momentum";

          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                  5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));

        end PartialCondition;

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={47,107,251})}));

    end InertDalton;

    package Translational
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector</html>"
      extends Modelica.Icons.Package;
      model Translational
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector, with efforts by default</html>"
        import FCSys.BaseClasses.Utilities.countTrue;
        import FCSys.BaseClasses.Utilities.enumerate;
        import FCSys.BaseClasses.Utilities.index;
        extends FCSys.BaseClasses.Icons.Conditions.Single;

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Axes with translational momentum included",
            compact=true));

        // Conditions
        replaceable Component.Velocity translationalX(final axis) if inclTransX
          constrainedby
          Conditions.ByConnector.Translational.Component.BaseClasses.PartialCondition(
            axis=Axis.x) "X-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransX),
          Placement(transformation(extent={{-52,32},{-32,52}})));
        replaceable Component.Velocity translationalY(final axis) if inclTransY
          constrainedby
          Conditions.ByConnector.Translational.Component.BaseClasses.PartialCondition(
            axis=Axis.y) "Y-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransY),
          Placement(transformation(extent={{-24,18},{-4,38}})));
        replaceable Component.Velocity translationalZ(final axis) if inclTransZ
          constrainedby
          Conditions.ByConnector.Translational.Component.BaseClasses.PartialCondition(
            axis=Axis.z) "Z-axis translational" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Conditions",enable=inclTransZ),
          Placement(transformation(extent={{4,4},{24,24}})));

        Connectors.Translational translational(final n_trans=n_trans)
          "Connector for advection or diffusion of translational momentum"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        Connectors.RealInputBus u
          "Input bus for values of specified conditions" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealOutputBus y "Output bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer cartTrans[n_trans]=index({inclTransX,
            inclTransY,inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

      equation
        // X-axis translational
        connect(translationalX.translational, translational) annotation (Line(
            points={{-42,38},{-42,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalX, translationalX.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,42},{-53,42}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalX.y, y.translationalX) annotation (Line(
            points={{-31,42},{90,42},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Y-axis translational
        connect(translationalY.translational, translational) annotation (Line(
            points={{-14,24},{-14,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalY, translationalY.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,28},{-25,28}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalY.y, y.translationalY) annotation (Line(
            points={{-3,28},{90,28},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

        // Z-axis translational
        connect(translationalZ.translational, translational) annotation (Line(
            points={{14,10},{14,-28},{5.55112e-16,-28},{5.55112e-16,-40}},
            color={239,142,1},
            smooth=Smooth.None));

        connect(u.translationalZ, translationalZ.u) annotation (Line(
            points={{-110,5.55112e-16},{-110,0},{-88,0},{-88,14},{3,14}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-2,3},{-2,3}}));
        connect(translationalZ.y, y.translationalZ) annotation (Line(
            points={{25,14},{90,14},{90,5.55112e-16},{110,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{2,3},{2,3}}));

      end Translational;

      model TranslationalIsolated
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector, with zero flows by default</html>"

        extends Translational(
          redeclare replaceable
            Conditions.ByConnector.Translational.Component.Force translationalX(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Translational.Component.Force translationalY(
              final axis),
          redeclare replaceable
            Conditions.ByConnector.Translational.Component.Force translationalZ(
              final axis));
        // See note in ChemicalReactionNoFlow.
        annotation (defaultComponentName="translational");

      end TranslationalIsolated;

      package Component "Conditions for a component of translational momentum"
        extends Modelica.Icons.Package;
        model Velocity "Specify velocity (measure force)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"),
            final y(final unit="l.m/T2") = translational.mPhidot[transCart[axis]]);

        equation
          if n_trans > 0 then
            translational.phi[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Velocity;

        model Force "Specify force (measure velocity)"
          extends BaseClasses.PartialCondition(
            final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l.m/T2"),
            final y(final unit="l/T") = translational.phi[transCart[axis]]);

        equation
          if n_trans > 0 then
            translational.mPhidot[transCart[axis]] = u_final;
          end if;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end Force;

        model Custom "Custom"
          extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
              y=translational.mPhidot[transCart[axis]]);

          Real x=translational.phi[transCart[axis]]
            "Expression to which the condition is applied"
            annotation (Dialog(group="Specification"));

        equation
          if n_trans > 0 then
            x = u_final;
          end if;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="translational",
            Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>chemical.phi[transCart[axis]]</code> and/or <code>chemical.mPhidot[transCart[axis]]</code>.</p></html>"));
        end Custom;

        package BaseClasses "Base classes (generally not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialCondition
            "Partial model for a translational condition"
            extends FCSys.BaseClasses.Icons.Conditions.Single;

            parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

            parameter Boolean internal=true "Use internal specification"
              annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(group="Specification"));

            replaceable Modelica.Blocks.Sources.RealExpression source if
              internal constrainedby Modelica.Blocks.Interfaces.SO
              "Source of internal specification" annotation (
              __Dymola_choicesFromPackage=true,
              Dialog(group="Specification",enable=internal),
              Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-70,30})));
            Connectors.RealInput u if not internal
              "Value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-110,0})));

            Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                  group="Measurement"), Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={110,0}), iconTransformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={110,0})));
            Connectors.Translational translational(final n_trans=n_trans)
              "Connector for advection or diffusion of translational momentum"
              annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

            constant ConditionType conditionType "Type of condition";
            // Note:  This is included so that the type of condition is recorded with
            // the results.

          protected
            outer parameter Integer n_trans
              "Number of components of translational momentum";
            outer parameter Integer cartTrans[:]
              "Cartesian-axis indices of the components of translational momentum";
            outer parameter Integer transCart[Axis]
              "Translational-momentum-component indices of the Cartesian axes";

            Connectors.RealOutputInternal u_final
              "Final value of specified condition" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={-20,0})));

          equation
            for i in 1:n_trans loop
              if cartTrans[i] <> axis then
                translational.mPhidot[i] = 0
                  "Not force in the other directions";
              end if;
            end for;
            connect(source.y, u_final) annotation (Line(
                points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},

                color={0,0,127},
                smooth=Smooth.None));

            connect(u, u_final) annotation (Line(
                points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                    5.55112e-16},{-20,5.55112e-16}},
                color={0,0,127},
                smooth=Smooth.None));
            annotation (defaultComponentName="translational");
          end PartialCondition;

          type ConditionType = enumeration(
              Velocity "Specify velocity (measure force)",
              Custom "Custom") "Types of conditions";

        end BaseClasses;

      end Component;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255})}));

    end Translational;

    package ThermalDiffusion
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.ThermalDiffusion\">ThermalDiffusion</a> connector</html>"
      extends Modelica.Icons.Package;

      model Temperature "Specify temperature (measure heat flow rate)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Temperature,
          u(final unit="l2.m/(N.T2)", displayUnit="K"),
          final y(unit="l2.m/T3") = thermal.Qdot,
          source(y=298.15*U.K));

      equation
        thermal.T = u_final;

      end Temperature;

      model HeatRate "Specify heat flow rate (measure temperature)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.HeatRate,
          u(final unit="l2.m/T3"),
          final y(
            final unit="l2.m/(N.T2)",
            displayUnit="K") = thermal.T);

      equation
        thermal.Qdot = u_final;

      end HeatRate;

      model Custom "Custom"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Custom,
          y=thermal.Qdot,
          source(y=298.15*U.K));

        Real x=thermal.T "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (defaultComponentName="thermal", Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>thermal.T</code> and/or <code>thermal.Qdot</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (generally not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a thermal condition"

          extends FCSys.BaseClasses.Icons.Conditions.Single;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression source if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification",enable=internal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,30})));

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with
          // the results.

          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-110,0})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                group="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));

          Connectors.ThermalDiffusion thermal
            "Connector to transport material, momentum, and energy of a single species"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        protected
          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-20,0})));

        equation
          connect(u, u_final) annotation (Line(
              points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{
                  -20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(source.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Specify temperature (measure heat flow rate)",
            HeatRate "Specify heat flow rate (measure temperature)",
            Custom "Custom") "Types of conditions";

      end BaseClasses;
      annotation (Icon(graphics={Ellipse(
              extent={{-70,50},{50,-70}},
              lineColor={170,0,0},
              fillPattern=FillPattern.Solid,
              fillColor={221,23,47})}));

    end ThermalDiffusion;
    annotation (Documentation(info="<html>
  <p>This package contains models to impose conditions on each of the declarative connectors
  established in <a href=\"modelica://FCSys.Connectors\">FCSys.Connectors</a>.  The subpackages
  are named according to the corresponding connector.</p>
</html>"));

  end ByConnector;

  record Environment "Environmental properties for a model"
    // extends FCSys.BaseClasses.Icons.Names.Top3;

    // Store the values of the base constants and units.
    final constant U.Bases.Base baseUnits=U.base "Base constants and units";

    parameter Boolean analysis=true "Include optional variables for analysis"
      annotation (choices(__Dymola_checkBox=true));

    parameter Q.PressureAbsolute p(nominal=U.atm) = U.atm "Pressure";
    parameter Q.TemperatureAbsolute T(nominal=300*U.K) = 298.15*U.K
      "Temperature";
    parameter Q.NumberAbsolute RH(displayUnit="%") = 1 "Relative humidity";
    parameter Q.NumberAbsolute n_O2_dry(
      final max=1,
      displayUnit="%") = 0.208
      "<html>Dry gas O<sub>2</sub> fraction (<i>n</i><sub>O2 dry</sub>)</html>";
    // Value from http://en.wikipedia.org/wiki/Oxygen
    parameter Q.Acceleration a[Axis]={0,Modelica.Constants.g_n*U.m/U.s^2,0}
      "Acceleration due to body forces";
    // The gravity component is positive because it's added to the transient
    // term in the Species model.
    parameter Q.PotentialLineic E[Axis]={0,0,0} "Electric field";
    final parameter Q.NumberAbsolute n_H2O(
      final max=1,
      displayUnit="%") = 0.2
      "<html>Gas H<sub>2</sub>O fraction (<i>n</i><sub>H2O</sub>)</html>";
    // TODO:  Cast this in terms of relative humidity.
    annotation (
      defaultComponentName="environment",
      defaultComponentPrefixes="inner",
      missingInnerMessage="
Your model is using an outer \"environment\" record, but an inner \"environment\"
record is not defined.  For simulation, drag FCSys.Conditions.Environment into
your model to specify global conditions and defaults.  Otherwise the default
settings will be used.
",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={
          Text(
            extent={{-120,60},{120,100}},
            textString="%name",
            lineColor={0,0,0}),
          Rectangle(
            extent={{-80,60},{80,-100}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-70,50},{70,-98}},
            lineColor={255,255,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={170,170,255}),
          Rectangle(
            extent={{-72,-60},{72,-100}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Line(points={{-70,-60},{70,-60}}, color={0,0,0}),
          Line(points={{-40,-20},{-10,-50},{40,0}}, color={0,0,0}),
          Ellipse(
            extent={{32,8},{48,-8}},
            pattern=LinePattern.None,
            lineColor={255,255,255},
            fillColor={50,50,50},
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
            lineColor={0,0,0}),
          Rectangle(
            extent={{-80,60},{80,-100}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash)}));

  end Environment;

  model Router "Connect two pairs of faces to pass through or cross over"
    extends FCSys.BaseClasses.Icons.Names.Top3;
    parameter Boolean crossOver=false "Cross over (otherwise, pass through)"
      annotation (choices(__Dymola_checkBox=true));
    Connectors.FaceBus negative1 "Negative face 1" annotation (Placement(
          transformation(extent={{-90,-50},{-70,-30}}, rotation=0),
          iconTransformation(extent={{-90,-50},{-70,-30}})));
    Connectors.FaceBus positive1 "Positive face 1" annotation (Placement(
          transformation(extent={{70,-50},{90,-30}}, rotation=0),
          iconTransformation(extent={{70,-50},{90,-30}})));
    Connectors.FaceBus negative2 "Negative face 2" annotation (Placement(
          transformation(extent={{-90,30},{-70,50}}, rotation=0),
          iconTransformation(extent={{-90,30},{-70,50}})));
    Connectors.FaceBus positive2 "Positive face 2" annotation (Placement(
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
          <img src=\"modelica://FCSys/Resources/Documentation/Conditions/Router/PassThrough.png\">
<br><b>a:</b>  Pass-through
        </td>
        <td align=center>
          <img src=\"modelica://FCSys/Resources/Documentation/Conditions/Router/CrossOver.png\">
<br><b>b:</b>  Cross-over
        </td>
      </tr>
      <tr align=center>
        <td colspan=2 align=center>Figure 1: Modes of connection.</td>
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

end Conditions;
