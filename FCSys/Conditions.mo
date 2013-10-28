within FCSys;
package Conditions "Models to specify and measure operating conditions"
  extends Modelica.Icons.SourcesPackage;
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;
    extends Modelica.Icons.UnderConstruction;

    model FaceCondition "Test the conditions for the face of a subregion"
      extends Modelica.Icons.Example;

      ByConnector.FaceBus.Single.Efforts face
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));
      Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclTransX=false,
        inclTransY=true,
        inclTransZ=false,
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
              "Resources/Scripts/Dymola/Conditions.Examples.FaceCondition.mos"
            "Conditions.Examples.FaceCondition.mos"));
    end FaceCondition;

    model FaceConditionPhases
      "Test the conditions for the face of a subregion with phases"
      import FCSys.Utilities.Coordinates.cartWrap;
      import Modelica.Math.BooleanVectors.countTrue;
      import Modelica.Math.BooleanVectors.enumerate;
      import Modelica.Math.BooleanVectors.index;
      extends Modelica.Icons.Example;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) =
        ones(3)*U.cm "Length" annotation (Dialog(group="Geometry",
            __Dymola_label="<html><b><i>L</i></b></html>"));
      final inner parameter Q.Volume V=product(L) "Volume";

      // Included faces
      parameter Boolean inclTransX=false "X" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));
      parameter Boolean inclTransY=true "Y" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));

      ByConnector.FaceBus.Single.Phases.Gas face(inclH2O=true, H2O(redeclare
            Conditions.ByConnector.Face.Single.Thermal.HeatRate thermal,
            redeclare Conditions.ByConnector.Face.Single.Material.Current
            material(source(y=U.A))))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));

      FCSys.Conditions.ByConnector.Amagat.Volume2 volume(n_phases=1)
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      FCSys.Phases.Gas gas(
        inclH2=false,
        inclH2O=true,
        final n_trans=n_trans,
        T_IC=environment.T,
        reduceThermal=false,
        H2O(T_IC=298.15*U.K))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    protected
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ} "true, if each pairs of faces is included";
      final inner parameter Boolean inclRot[Axis]={inclTransY and inclTransZ,
          inclTransZ and inclTransX,inclTransX and inclTransY}
        "true, if each axis of rotation has all its tangential faces included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartTrans[n_trans]=index(inclTrans)
        "Cartesian-axis indices of the transport axes";
      final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
        "Transport-axis indices of the Cartesian axes";

    equation
      connect(gas.yPositive, face.face) annotation (Line(
          points={{0,10},{0,16},{6.10623e-016,16},{6.10623e-016,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(volume.dalton[1], gas.dalton) annotation (Line(
          points={{11,-11},{8,-8}},
          color={47,107,251},
          smooth=Smooth.None));
      annotation (experiment);
    end FaceConditionPhases;

    model Router
      "<html>Test the <a href=\"modelica://FCSys.Conditions.Router\">Router</a> model</html>"
      extends Modelica.Icons.Example;

      Conditions.Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      ByConnector.FaceBus.Single.Flows fastFlow(gas(inclH2=true, H2(materialSet(
                y=-U.A)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-30,20})));
      ByConnector.FaceBus.Single.Flows slowFlow(gas(inclH2=true, H2(materialSet(
                y=-U.mA)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-30,-20})));

      ByConnector.FaceBus.Single.Efforts sink2(gas(inclH2=true)) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-20})));
      ByConnector.FaceBus.Single.Efforts sink1(gas(inclH2=true)) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,20})));
    equation
      connect(router.positive2, sink1.face) annotation (Line(
          points={{8,4},{20,4},{20,20},{26,20},{26,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(router.positive1, sink2.face) annotation (Line(
          points={{8,-4},{20,-4},{20,-20},{26,-20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(fastFlow.face, router.negative2) annotation (Line(
          points={{-26,20},{-20,20},{-20,4},{-8,4}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(slowFlow.face, router.negative1) annotation (Line(
          points={{-26,-20},{-20,-20},{-20,-4},{-8,-4}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics), Commands(file=
              "Resources/Scripts/Dymola/Conditions.Examples.Router.mos"
            "Conditions.Examples.Router.mos"));
    end Router;

    model AnodeAdapter
      "<html>Test the <a href=\"modelica://FCSys.Conditions.Adapters.Anode\">Anode</a> adapter</html>"
      extends Modelica.Icons.Example;

      inner Modelica.Fluid.System system(T_ambient=293.15 + 5)
        annotation (Placement(transformation(extent={{40,70},{60,90}})));
      inner Conditions.Environment environment(T=350*U.K)
        annotation (Placement(transformation(extent={{70,72},{90,92}})));
      FCSys.Subregions.SubregionNoIonomer subregion(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        gas(inclH2=true, inclH2O=true),
        graphite('inclC+'=true, 'incle-'=true),
        liquid(inclH2O=true))
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Adapters.MSL.Anode anodeAdapter(redeclare package LiquidMedium =
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
              "Resources/Scripts/Dymola/Conditions.Examples.Adapteminus.mos"
            "Conditions.Examples.Adapteminus.mos"));
    end AnodeAdapter;

    model Stoichiometry "Test the stoichiometry of a reaction"

      extends Modelica.Icons.Example;
      import Modelica.Math.BooleanVectors.countTrue;

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));
      parameter Boolean inclTransY=false "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Included transport axes",
          compact=true));

      inner Conditions.Environment environment(T=360*U.K)
        annotation (Placement(transformation(extent={{40,40},{60,60}})));

      replaceable Conditions.ByConnector.Electrochemical.Current speciesA(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=2000*U.K,
        redeclare Modelica.Blocks.Sources.Ramp source(duration=100, height=-1*U.A))
        annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
      replaceable Conditions.ByConnector.Electrochemical.Potential speciesB(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=3000*U.K)
        annotation (Placement(transformation(extent={{30,0},{50,20}})));
      replaceable Conditions.ByConnector.Electrochemical.Potential speciesC(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=4000*U.K)
        annotation (Placement(transformation(extent={{30,-20},{50,0}})));
      FCSys.Conditions.Adapters.ElectrochemicalReaction A(
        m=U.g/U.mol,
        final n_trans=n_trans,
        n=-1)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      FCSys.Conditions.Adapters.ElectrochemicalReaction B(
        m=U.g/U.mol,
        final n_trans=n_trans,
        n=2) annotation (Placement(transformation(extent={{30,-10},{10,10}})));

      FCSys.Conditions.Adapters.ElectrochemicalReaction C(
        m=U.g/U.mol,
        final n_trans=n_trans,
        n=2) annotation (Placement(transformation(extent={{30,-30},{10,-10}})));

    protected
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";

    equation
      connect(A.electrochemical, speciesA.electrochemical) annotation (Line(
          points={{-24,0},{-40,0},{-40,6}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(B.electrochemical, speciesB.electrochemical) annotation (Line(
          points={{24,0},{40,0},{40,6}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(C.electrochemical, speciesC.electrochemical) annotation (Line(
          points={{24,-20},{40,-20},{40,-14}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(C.reaction, B.reaction) annotation (Line(
          points={{16,-20},{0,-20},{0,0},{16,0}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(B.reaction, A.reaction) annotation (Line(
          points={{16,0},{-16,0}},
          color={221,23,47},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=120),
        Commands(file=
              "Resources/Scripts/Dymola/Conditions.Examples.Stoichiometry.mos"
            "Conditions.Examples.Stoichiometry.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        __Dymola_experimentSetupOutput);
    end Stoichiometry;
  end Examples;

  package Adapters
    "<html>Interfaces to the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
    extends Modelica.Icons.Package;
    extends Modelica.Icons.UnderConstruction;

    model AmagatDalton
      "<html>Adapter between the <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> and <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connectors</html>"

      extends FCSys.Icons.Names.Top1;

      FCSys.Connectors.Amagat amagat "Connector for additivity of volume"
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}}),
            iconTransformation(extent={{-50,-10},{-30,10}})));
      FCSys.Connectors.Dalton dalton "Connector for additivity of pressure"
        annotation (Placement(transformation(extent={{10,-10},{30,10}}),
            iconTransformation(extent={{30,-10},{50,10}})));

    equation
      // Static balances
      0 = amagat.p + dalton.p "Pressure";
      0 = amagat.V + dalton.V "Volume";

      annotation (
        Documentation(info="<html><p>This model is used to convert the representation of mixtures 
    between Amagat's law of partial volumes and Dalton's law of partial pressures.</p>
    
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                100,100}}), graphics={Line(
                  points={{-30,0},{30,0}},
                  color={47,107,251},
                  smooth=Smooth.None),Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={47,107,251},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}));
    end AmagatDalton;

    model ElectrochemicalReaction
      "<html>Adapter between the <a href=\"modelica://FCSys.Connectors.Electrochemical\">Electrochemical</a> and <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connectors</html>"

      //extends FCSys.Icons.Names.Top1;

      constant Integer n_trans(min=1,max=3)
        "Number of components of translational momentum" annotation (Dialog(
            __Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
      // Note:  This must be a constant rather than a parameter due to errors
      // in Dymola 2014.
      parameter Integer n "Stoichiometric coefficient"
        annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
      parameter Q.MassSpecific m "Specific mass" annotation (Dialog(group=
              "Material properties", __Dymola_label="<html><i>m</i></html>"));

      // Auxiliary variables (for analysis)
      /*
  output Q.Velocity phi[n_trans](each stateSelect=StateSelect.never) = 
    actualStream(electrochemical.phi) if environment.analysis "Velocity of the stream";
  output Q.PotentialAbsolute sT(stateSelect=StateSelect.never) = actualStream(
    electrochemical.sT) if environment.analysis 
    "Specific entropy-temperature product of the stream";
*/

      Connectors.Electrochemical electrochemical(redeclare final constant
          Integer n_trans=n_trans)
        "Connector for a species in a chemical reaction" annotation (Placement(
            transformation(extent={{-30,-10},{-10,10}}), iconTransformation(
              extent={{-50,-10},{-30,10}})));
      // Note:  This redeclaration is necessary due to errors in Dymola 2014.
      Connectors.Reaction reaction(final n_trans=n_trans)
        "Connector for an electrochemical reaction" annotation (Placement(
            transformation(extent={{10,-10},{30,10}}), iconTransformation(
              extent={{30,-10},{50,10}})));

    equation
      // Equal intensive properties
      reaction.w = n*electrochemical.w "Electrochemical potential";
      reaction.phi = electrochemical.phi "Velocity (upon outflow)";
      reaction.sT = electrochemical.sT
        "Specific entropy-temperature product (upon outflow)";

      // Conservation (without storage)
      0 = electrochemical.Ndot + n*reaction.Ndot "Material";
      zeros(n_trans) = m*actualStream(electrochemical.phi)*electrochemical.Ndot
         + reaction.mPhidot "Translational momentum";
      0 = actualStream(electrochemical.sT)*electrochemical.Ndot + reaction.Qdot
        "Energy";
      annotation (
        Documentation(info="<html><p>This model is used to add the stoichiometrically-weighted electrochemical potential
    of a species to the net electrochemical potential of a reaction.  The species is produced at the
    stoichiometrically-weighted rate of the reaction.</p>
    
    <p>
    For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={Line(
                  points={{-30,0},{30,0}},
                  color={221,23,47},
                  smooth=Smooth.None),Text(
                  extent={{-100,20},{100,60}},
                  lineColor={0,0,0},
                  textString="%n %name"),Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={221,23,47},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end ElectrochemicalReaction;

    package MSL
      "<html>Adapters to the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
      extends Modelica.Icons.Package;

      // TODO: Create a wrapper for the whole cell.
      model Anode
        "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the face connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>"
        extends FCSys.Icons.Names.Top4;

        replaceable package GasMedium = Media.AnodeGas constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the gas"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));
        replaceable package LiquidMedium =
            Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Connectors.FaceBus face
          "Multi-species connector for translational momentum and heat"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
            Medium = GasMedium) "Modelica fluid port for the gas" annotation (
            Placement(transformation(extent={{70,70},{90,90}}),
              iconTransformation(extent={{70,50},{90,70}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,10},{90,30}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -14},{90,6}}), iconTransformation(extent={{70,-30},{90,-10}})));
        Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final
            package Medium = LiquidMedium) "Modelica fluid port for the liquid"
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
            points={{8,6.10623e-16},{24,6.10623e-16},{24,0},{40,0},{40,-80},{80,
                -80}},
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
        extends FCSys.Icons.Names.Top4;

        replaceable package GasMedium = Adapters.Media.CathodeGas
          constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model for the gas" annotation (choicesAllMatching=true,
            Dialog(group="Material properties"));
        replaceable package LiquidMedium =
            Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Connectors.FaceBus face
          "Multi-species connector for translational momentum and heat"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
            Medium = GasMedium) "Modelica fluid port for the gas" annotation (
            Placement(transformation(extent={{70,70},{90,90}}),
              iconTransformation(extent={{70,50},{90,70}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,10},{90,30}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -14},{90,6}}), iconTransformation(extent={{70,-30},{90,-10}})));
        Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final
            package Medium = LiquidMedium) "Modelica fluid port for the liquid"
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
            points={{8,6.10623e-16},{24,6.10623e-16},{24,0},{40,0},{40,-80},{80,
                -80}},
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

      model Conductor
        "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the face connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>, with only the graphite phase included"
        extends FCSys.Icons.Names.Top4;

        replaceable package GasMedium = Adapters.Media.CathodeGas
          constrainedby Modelica.Media.Interfaces.PartialMedium
          "Medium model for the gas" annotation (choicesAllMatching=true,
            Dialog(group="Material properties"));
        replaceable package LiquidMedium =
            Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Connectors.FaceBus face
          "Multi-species connector for translational momentum and heat"
          annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{70,
                  -50},{90,-30}}), iconTransformation(extent={{70,-50},{90,-30}})));

        Phases.Graphite graphite "Graphite subadapter"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
          "Modelica translational flanges" annotation (Placement(transformation(
                extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));

      equation
        connect(graphite.face, face.graphite) annotation (Line(
            points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5));

        connect(graphite.pin, pin) annotation (Line(
            points={{8,4},{40,4},{40,40},{80,40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(graphite.heatPort, heatPort) annotation (Line(
            points={{8,-4},{40,-4},{40,-40},{80,-40}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(graphite.flange, flange) annotation (Line(
            points={{8,6.10623e-16},{24,6.10623e-16},{24,0},{40,0},{40,
                5.55112e-16},{80,5.55112e-16}},
            color={0,127,0},
            smooth=Smooth.None));

        annotation (Icon(graphics={Line(
                      points={{0,40},{0,-40}},
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
                      points={{0,-40},{80,-40}},
                      color={191,0,0},
                      smooth=Smooth.None),Line(
                      points={{0,0},{70,0}},
                      color={0,127,0},
                      smooth=Smooth.None)}), Diagram(graphics));
      end Conductor;

      package Phases "Adapters for material phases"
        extends Modelica.Icons.Package;

        model AnodeGas
          "<html>Adapter for PEMFC anode gas between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends Partial;

          replaceable package Medium = Media.AnodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Species.FluidNeutral H2(redeclare package Medium =
                Modelica.Media.IdealGases.SingleGases.H2 (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false), redeclare package Data =
                Characteristics.H2.Gas)
            annotation (Placement(transformation(extent={{-10,10},{10,30}})));
          Species.FluidNeutral H2O(redeclare package Data =
                Characteristics.H2O.Gas (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false), redeclare final package
              Medium = Modelica.Media.IdealGases.SingleGases.H2O)
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
          Junctions.Junction2 junction
            annotation (Placement(transformation(extent={{62,30},{42,50}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));

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
          extends Partial;

          replaceable package Medium = Media.CathodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Junctions.Junction3 junction(
            redeclare package Medium1 =
                Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false),
            redeclare package Medium2 =
                Modelica.Media.IdealGases.SingleGases.N2 (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false),
            redeclare package Medium3 =
                Modelica.Media.IdealGases.SingleGases.O2 (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false),
            redeclare package MixtureMedium = Medium)
            annotation (Placement(transformation(extent={{60,30},{40,50}})));

          Conditions.Adapters.MSL.Species.FluidNeutral H2O(redeclare package
              Data = Characteristics.H2O.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false))
            annotation (Placement(transformation(extent={{-10,10},{10,30}})));
          Conditions.Adapters.MSL.Species.FluidNeutral N2(redeclare package
              Data = Characteristics.N2.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.N2 (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false))
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
          Conditions.Adapters.MSL.Species.FluidNeutral O2(redeclare package
              Data = Characteristics.O2.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.O2 (referenceChoice=
                    Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                  excludeEnthalpyOfFormation=false))
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));

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
          extends Partial;

          Species.'e-' 'e-'(redeclare package Data =
                Characteristics.'e-'.Graphite)
            annotation (Placement(transformation(extent={{-10,30},{10,50}})));

          Species.Solid 'C+'(redeclare package Data =
                Characteristics.'C+'.Graphite)
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
          Modelica.Electrical.Analog.Interfaces.NegativePin pin
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{70,30},{90,50}}), iconTransformation(extent={{70,30},
                    {90,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));

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
          extends Partial;

          replaceable package Medium =
              Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "Medium model (Modelica)" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Conditions.Adapters.MSL.Species.FluidNeutral H2O(redeclare package
              Data = Characteristics.H2O.Liquid, redeclare final package Medium
              = Medium)
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));

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

      protected
        partial model Partial
          "<html>Base model for adapter for a phase between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.Icons.Names.Top3;

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

        end Partial;
      end Phases;

      package Species "Adapters for single species"
        extends Modelica.Icons.Package;

        model 'e-'
          "<html>Adapter to connect e<sup>-</sup> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (electrical and heat only)</html>"
          import FCSys.Utilities.inSign;
          extends FCSys.Icons.Names.Top2;

          // Geometry
          parameter Q.Area A=U.cm^2 "Area of the interface" annotation (Dialog(
                group="Geometry", __Dymola_label="<html><i>A</i></html>"));
          parameter Side side=Side.n "FCSys side of the interface"
            annotation (Dialog(group="Geometry"));

          replaceable package Data = Characteristics.'e-'.Graphite
            constrainedby Characteristics.BaseClasses.Characteristic
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
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{70,-10},{90,10}}), iconTransformation(extent={{70,-10},
                    {90,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-50},{90,-30}}), iconTransformation(extent={{70,-50},{90,
                    -30}})));

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
          import FCSys.Utilities.inSign;
          extends FCSys.Icons.Names.Top3;

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
          replaceable package Medium =
              Modelica.Media.IdealGases.SingleGases.H2O constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "Medium model (for Modelica)" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Medium.BaseProperties medium "Base properties of the fluid";
          Q.Current I "Material current";

          Connectors.Face face
            "Connector for material, momentum, and energy of a single species"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));
          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-50},{90,-30}}), iconTransformation(extent={{70,-50},{90,
                    -30}})));
          Modelica.Electrical.Analog.Interfaces.NegativePin pin
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{70,-90},{90,-70}}), iconTransformation(extent={{70,-90},
                    {90,-70}})));

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
          import assert = FCSys.Utilities.assertEval;
          extends FCSys.Icons.Names.Top3;

          parameter Q.Area A=U.cm^2 "Area of the interface" annotation (Dialog(
                group="Geometry", __Dymola_label="<html><i>A</i></html>"));
          replaceable package Data = Characteristics.BaseClasses.Characteristic
            "Characteristic data (for FCSys)" annotation (
            Dialog(group="Material properties"),
            __Dymola_choicesAllMatching=true,
            Placement(transformation(extent={{-60,40},{-40,60}}),
                iconTransformation(extent={{-10,90},{10,110}})));
          replaceable package Medium =
              Modelica.Media.IdealGases.SingleGases.H2O constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "Medium model (for Modelica)" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Medium.BaseProperties medium "Base properties of the fluid";

          Connectors.Face face
            "Connector for material, momentum, and energy of a single species"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));
          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Axis]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{70,-10},{90,10}}), iconTransformation(
                  extent={{70,-10},{90,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-50},{90,-30}}), iconTransformation(extent={{70,-50},{90,
                    -30}})));

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
          0 = face.Ndot + face.phi[1]*face.rho*A + (fluidPort.m_flow/medium.MM)
            *U.mol/U.s "Material";
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
          extends FCSys.Icons.Names.Top2;

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
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-50},{90,-30}}), iconTransformation(extent={{70,-50},{90,
                    -30}})));

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

        model Junction2
          "Junction between two pure substances and their mixture"
          import assert = FCSys.Utilities.assertEval;
          extends Partial;

          replaceable package Medium1 =
              Modelica.Media.IdealGases.SingleGases.H2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false) constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium2 =
              Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false) constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final
              package Medium = Medium1) "Fluid port for the 1st pure substance"
            annotation (Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final
              package Medium = Medium2) "Fluid port for the 2nd pure substance"
            annotation (Placement(transformation(extent={{70,-50},{90,-30}}),
                iconTransformation(extent={{70,-50},{90,-30}})));

        initial equation
          // Check the number and names of substances
          assert(MixtureMedium.nS == 2,
            "The mixture medium must have exactly two substances.");
          assert(MixtureMedium.substanceNames[1] == Medium1.substanceNames[1],
            "The first substance of the mixture medium (MixtureMedium) is \""
             + MixtureMedium.substanceNames[1] + "\",
but the first pure substance is \"" + Medium1.substanceNames[1] + "\".");
          assert(MixtureMedium.substanceNames[2] == Medium2.substanceNames[1],
            "The second substance of the mixture medium (MixtureMedium) is \""
             + MixtureMedium.substanceNames[2] + "\",
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
          import assert = FCSys.Utilities.assertEval;
          extends Partial(redeclare replaceable package MixtureMedium =
                Media.CathodeGas);

          replaceable package Medium1 =
              Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false) constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium2 =
              Modelica.Media.IdealGases.SingleGases.N2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false) constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium3 =
              Modelica.Media.IdealGases.SingleGases.O2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false) constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 3<sup>rd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final
              package Medium = Medium1) "Fluid port for the 1st pure substance"
            annotation (Placement(transformation(extent={{70,30},{90,50}}),
                iconTransformation(extent={{70,30},{90,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final
              package Medium = Medium2) "Fluid port for the 2nd pure substance"
            annotation (Placement(transformation(extent={{70,-10},{90,10}}),
                iconTransformation(extent={{70,-10},{90,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort3(redeclare final
              package Medium = Medium3) "Fluid port for the 3rd pure substance"
            annotation (Placement(transformation(extent={{70,-50},{90,-30}}),
                iconTransformation(extent={{70,-50},{90,-30}})));

        initial equation
          // Check the number and names of substances
          assert(MixtureMedium.nS == 3,
            "The mixture medium must have exactly three substances.");
          assert(MixtureMedium.substanceNames[1] == Medium1.substanceNames[1],
            "The first substance of the mixture medium (MixtureMedium) is \""
             + MixtureMedium.substanceNames[1] + "\",
but the first pure substance is \"" + Medium1.substanceNames[1] + "\".");
          assert(MixtureMedium.substanceNames[2] == Medium2.substanceNames[1],
            "The second substance of the mixture medium (MixtureMedium) is \""
             + MixtureMedium.substanceNames[2] + "\",
but the second pure substance is \"" + Medium2.substanceNames[1] + "\".");
          assert(MixtureMedium.substanceNames[3] == Medium3.substanceNames[1],
            "The second substance of the mixture medium (MixtureMedium) is \""
             + MixtureMedium.substanceNames[3] + "\",
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

      protected
        partial model Partial
          "Base model for a junction between pure substances and their mixture"
          extends FCSys.Icons.Names.Top3;

          replaceable package MixtureMedium = Media.AnodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium
            "Medium model for the mixture" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_a mixturePort(redeclare final
              package Medium = MixtureMedium) "Fluid port for the mixture"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));

          FCSys.Conditions.Adapters.MSL.SIunits.MassFraction X[MixtureMedium.nX]
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
        end Partial;
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
    end MSL;

  end Adapters;

  package ByConnector "Conditions for each type of connector"
    package Electrostatic
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Electrostatic\">Electrostatic</a> connector</html>"
      extends Modelica.Icons.Package;

      model Potential "Specify potential (measure charge)"
        extends Partial(final y=electrostatic.Z);

      equation
        electrostatic.w = u_final;

      end Potential;

      model Charge "Specify charge (measure potential)"
        extends Partial(final y=electrostatic.w);

      equation
        electrostatic.Z = u_final;

      end Charge;

      partial model Partial "Base model for a electrostatic condition"

        extends FCSys.Icons.Conditions.SingleShort;

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

        Connectors.RealOutput y "Measurement expression" annotation (Dialog(tab
              ="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.Electrostatic electrostatic "Static electrical connector"
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
            points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{-20,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (
          defaultComponentName="electrostatic",
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end Partial;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={239,142,1},
              fillPattern=FillPattern.Solid,
              fillColor={255,195,38})}));
    end Electrostatic;
    extends Modelica.Icons.Package;

    package Electrochemical
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Electrochemical\">Electrochemical</a> connector</html>"
      extends Modelica.Icons.Package;

      model Potential "Specify electrochemical potential (measure current)"
        extends Partial(final y=electrochemical.Ndot);

      equation
        electrochemical.w = u_final;

      end Potential;

      model Current "Specify current (measure electrochemical potential)"
        extends Partial(final y=electrochemical.w);

      equation
        electrochemical.Ndot = u_final;

        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent=
                  {{-100,-100},{100,100}}), graphics));
      end Current;

      partial model Partial "Base model for a material condition"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.index;
        extends FCSys.Icons.Conditions.SingleShort;

        parameter Boolean internal=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Specification of material condition"));

        replaceable Modelica.Blocks.Sources.RealExpression source if internal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Specification of material condition", enable=internal),

          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,10})));

        // Properties upon outflow
        parameter Q.Velocity phi[Axis]={0,0,0} "Velocity" annotation (Dialog(
              group="Properties upon outflow", __Dymola_label=
                "<html><i><b>&phi;</b></i></html>"));
        parameter Q.PotentialAbsolute sT=3000*U.K
          "Specific entropy-temperature product" annotation (Dialog(group=
                "Properties upon outflow", __Dymola_label=
                "<html><i>sT</i></html>"));

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));

        Connectors.RealInput u if not internal "Value of specified condition"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));

        Connectors.RealOutput y "Measurement expression" annotation (Dialog(tab
              ="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));
        output Q.Velocity phi_actual[n_trans]=actualStream(electrochemical.phi)
          "Velocity of the actual stream";
        output Q.Potential sT_actual=actualStream(electrochemical.sT)
          "Specific entropy-temperature product of the actual stream";
        Connectors.Electrochemical electrochemical(final n_trans=n_trans)
          "Connector for a species of a chemical reaction"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        // final formula=formula

      protected
        final parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final parameter Integer cartTrans[n_trans]=index({inclTransX,inclTransY,
            inclTransZ})
          "Cartesian-axis indices of the components of translational momentum";
        Connectors.RealOutputInternal u_final
          "Final value of specified condition" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,0}),iconTransformation(extent={{-10,-10},{10,10}},
                origin={-20,0})));

      equation
        electrochemical.phi = phi[cartTrans];
        electrochemical.sT = sT;

        connect(source.y, u_final) annotation (Line(
            points={{-69,10},{-60,10},{-60,5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(u, u_final) annotation (Line(
            points={{-110,5.55112e-16},{-88,0},{-66,1.11022e-15},{-66,
                5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (defaultComponentName="chemical");
      end Partial;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={170,0,0},
              fillPattern=FillPattern.Solid,
              fillColor={221,23,47})}));

    end Electrochemical;

    package Reaction
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector</html>"
      extends Modelica.Icons.Package;

      model Flows
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.ElectrochemNegative\">ElectrochemNegative</a> or <a href=\"modelica://FCSys.Connectors.ElectrochemPositive\">ElectrochemPositive</a> connector, with flows specified by default</html>"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.enumerate;
        import Modelica.Blocks.Sources;
        extends FCSys.Icons.Conditions.Single;

        // Specification
        // -------------
        // Material
        replaceable function materialSpec = Material.potential constrainedby
          Material.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          choicesAllMatching=true,
          Dialog(tab="Specification", group="Material"));
        parameter Boolean internalMaterial=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(tab="Specification", group="Material"));
        replaceable Sources.RealExpression materialSet if internalMaterial
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Material",
            enable=internalMaterial),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,90})));

        //
        // X-axis translational
        replaceable function transXSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX),
          Placement(transformation(extent={{-52,18},{-32,38}})));

        parameter Boolean internalTransX=true if inclTransX
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX));
        replaceable Sources.RealExpression transXSet if inclTransX and
          internalTransX constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX and internalTransX),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,50})));

        //
        // Y-axis translational
        replaceable function transYSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY),
          Placement(transformation(extent={{-24,4},{-4,24}})));

        parameter Boolean internalTransY=true if inclTransY
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY));
        replaceable Sources.RealExpression transYSet if inclTransY and
          internalTransY constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY and internalTransY),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,10})));

        //
        // Z-axis translational
        replaceable function transZSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalTransZ=true if inclTransZ
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ));
        replaceable Sources.RealExpression transZSet if inclTransZ and
          internalTransZ constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ and internalTransZ),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-30})));

        //
        // Thermal
        replaceable function thermalSpec = Thermal.heatRate constrainedby
          Thermal.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Specification", group="Thermal"),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalThermal=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(tab="Specification", group="Thermal"));
        replaceable Sources.RealExpression thermalSet if internalThermal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Thermal",
            enable=internalThermal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-70})));

        // Measurement
        // -----------
        // Material
        replaceable function materialMeas = Material.reactionRate
          constrainedby Material.Partial "Material quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(group="Measurement"));

        // X-axis translational
        replaceable function transXMeas = Translational.velocity constrainedby
          Translational.Partial "X-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(group="Measurement"));

        // Y-axis translational
        replaceable function transYMeas = Translational.velocity constrainedby
          Translational.Partial "Y-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(group="Measurement"));

        // Z-axis translational
        replaceable function transZMeas = Translational.velocity constrainedby
          Translational.Partial "Z-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(group="Measurement"));

        // Thermal
        replaceable function thermalMeas = Thermal.specificEntropyTemperature
          constrainedby Thermal.Partial "Thermal quantity" annotation (
            __Dymola_choicesFromPackage=true,Dialog(tab="Measurement"));

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));

        // Inputs
        Connectors.RealInput u_material if not internalMaterial
          "Material specification" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,80})));
        Connectors.RealInput u_transX if inclTransX and not internalTransX
          "X-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,40})));
        Connectors.RealInput u_transY if inclTransY and not internalTransY
          "Y-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealInput u_transZ if inclTransZ and not internalTransZ
          "Z-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-40})));
        Connectors.RealInput u_thermal if not internalThermal
          "Thermal specification" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-80})));

        // Outputs
        final Connectors.RealOutput y_material=materialMeas(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot) "Material measurement" annotation (Dialog(tab=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,80}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,80})));
        final Connectors.RealOutput y_transX=transXMeas(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));
        final Connectors.RealOutput y_transY=transYMeas(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));
        final Connectors.RealOutput y_transZ=transZMeas(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40})));
        final Connectors.RealOutput y_thermal=thermalMeas(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot) "Thermal measurement" annotation (Dialog(group
              ="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-80}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-80})));

        Connectors.Reaction reaction(final n_trans=n_trans)
          "Stoichiometric connector for the electrochemical reaction"
          annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

        Connectors.RealOutputInternal _u_material=materialSpec(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot)
          "Internal, working value of material specification" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,80})));
        Connectors.RealOutputInternal _u_transX=transXSpec(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "Internal, working value of X-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,40})));
        Connectors.RealOutputInternal _u_transY=transYSpec(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Internal, working value of Y-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,0})));
        Connectors.RealOutputInternal _u_transZ=transZSpec(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Internal, working value of Z-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-40})));
        Connectors.RealOutputInternal _u_thermal=thermalSpec(
                  reaction.Ndot,
                  reaction.w,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot)
          "Internal, working value of thermal specification" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-80})));

      equation
        // Material
        connect(u_material, _u_material) annotation (Line(
            points={{-110,80},{-36,80}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(materialSet.y, _u_material) annotation (Line(
            points={{-69,90},{-60,90},{-60,80},{-36,80}},
            color={0,0,127},
            smooth=Smooth.None));

        // X-axis translational
        connect(u_transX, _u_transX) annotation (Line(
            points={{-110,40},{-36,40}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transXSet.y, _u_transX) annotation (Line(
            points={{-69,50},{-60,50},{-60,40},{-36,40}},
            color={0,0,127},
            smooth=Smooth.None));

        // Y-axis translational
        connect(u_transY, _u_transY) annotation (Line(
            points={{-110,5.55112e-16},{-36,5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transYSet.y, _u_transY) annotation (Line(
            points={{-69,10},{-60,10},{-60,5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        // Z-axis translational
        connect(u_transZ, _u_transZ) annotation (Line(
            points={{-110,-40},{-36,-40}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transZSet.y, _u_transZ) annotation (Line(
            points={{-69,-30},{-60,-30},{-60,-40},{-36,-40}},
            color={0,0,127},
            smooth=Smooth.None));

        // Thermal
        connect(u_thermal, _u_thermal) annotation (Line(
            points={{-110,-80},{-36,-80}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(thermalSet.y, _u_thermal) annotation (Line(
            points={{-69,-70},{-60,-70},{-60,-80},{-36,-80}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (
          defaultComponentName="reaction",
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics),
          Icon(graphics));
      end Flows;

      model Efforts
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.ElectrochemNegative\">ElectrochemNegative</a> or <a href=\"modelica://FCSys.Connectors.ElectrochemPositive\">ElectrochemPositive</a> connector, with efforts specified by default</html>"
        extends FCSys.Conditions.ByConnector.Reaction.Flows(
          redeclare replaceable function materialSpec =
              Conditions.ByConnector.Reaction.Material.reactionRate,
          redeclare replaceable function transXSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function transYSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function transZSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function thermalSpec =
              Conditions.ByConnector.Reaction.Thermal.specificEntropyTemperature,

          redeclare replaceable function materialMeas =
              Conditions.ByConnector.Reaction.Material.potential,
          redeclare replaceable function transXMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function transYMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function transZMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function thermalMeas =
              Conditions.ByConnector.Reaction.Thermal.heatRate);

        // Note:  In Dymola 2014, the paths must be explicitly given to prevent
        // the error "Cannot show the parameter dialog for redeclared class [...]".
        annotation (defaultComponentName="reaction");

      end Efforts;

      package Material "Material conditions"
        extends Modelica.Icons.Package;

        function reactionRate "Reaction rate"
          extends Partial;

        algorithm
          x := Ndot;
          annotation (Inline=true);
        end reactionRate;

        function potential "Electrochemical potential"
          extends Partial;

        algorithm
          x := g;
          annotation (Inline=true);
        end potential;

        partial function Partial
          "Template of a function to select a material quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "<html>Rate of reaction (<i>N&#775;</i>)</html>";
          input Q.Potential g
            "<html>Electrochemical potential (<i>g</i>)</html>";

          // Translational advection
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal advection
          input Q.PotentialAbsolute sT "Specific entropy-temperature product";
          input Q.Power Qdot
            "<html>Rate of thermal advection (<i>Q&#775;</i>)</html>";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Material;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;

        function velocity "Velocity"
          extends Partial;

        algorithm
          x := phi[i];
          annotation (Inline=true);
        end velocity;

        function force "Force"
          extends Partial;

        algorithm
          x := mPhidot[i];
          annotation (Inline=true);
        end force;

        partial function Partial
          "Template of a function to select a translational quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "<html>Rate of reaction (<i>N&#775;</i>)</html>";
          input Q.Potential g
            "<html>Electrochemical potential (<i>g</i>)</html>";

          // Translational advection
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal advection
          input Q.PotentialAbsolute sT "Specific entropy-temperature product";
          input Q.Power Qdot
            "<html>Rate of thermal advection (<i>Q&#775;</i>)</html>";

          input Integer i(min=1,max=3) "Index of the translational axis";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Translational;

      package Thermal "Thermal conditions"
        extends Modelica.Icons.Package;

        function specificEntropyTemperature
          "Specific entropy-temperature product"
          extends Partial;

        algorithm
          x := sT;
          annotation (Inline=true);
        end specificEntropyTemperature;

        function heatRate "Heat flow rate"
          extends Partial;

        algorithm
          x := Qdot;
          annotation (Inline=true);
        end heatRate;

        partial function Partial
          "Template of a function to select a thermal quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "<html>Rate of reaction (<i>N&#775;</i>)</html>";
          input Q.Potential g
            "<html>Electrochemical potential (<i>g</i>)</html>";

          // Translational advection
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal advection
          input Q.PotentialAbsolute sT "Specific entropy-temperature product";
          input Q.Power Qdot
            "<html>Rate of thermal advection (<i>Q&#775;</i>)</html>";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Thermal;

      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={170,0,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}), Ellipse(
              extent={{-30,30},{30,-30}},
              fillColor={221,23,47},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0})}));
    end Reaction;

    package FaceBus
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
      extends Modelica.Icons.Package;

      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"

        extends Modelica.Icons.Package;

        model Flows
          "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors, with flows specified by default</html>"

          extends FCSys.Icons.Conditions.SingleShort;

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
          annotation (defaultComponentName="face", Diagram(graphics));
        end Flows;

        model FlowsFluidOnly
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows specified by default and only the fluid phases included</html>"

          extends FCSys.Icons.Conditions.SingleShort;

          Phases.Gas gas "Gas" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Phases.Liquid liquid "Liquid" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-110,0})));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={110,0})));

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
              points={{5.55112e-16,50},{6.10623e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(gas.y, y.gas) annotation (Line(
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
              points={{5.55112e-16,50},{6.10623e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(liquid.y, y.liquid) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (defaultComponentName="face");
        end FlowsFluidOnly;

        model FlowsGraphiteOnly
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows specified by default and only the graphite phase</html>"

          extends FCSys.Icons.Conditions.SingleShort;

          Phases.Graphite graphite "Graphite" annotation (Dialog(group=
                  "Phases (click to edit)", __Dymola_descriptionLabel=true),
              Placement(transformation(extent={{-10,-10},{10,10}})));

          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50}), iconTransformation(extent={{-10,-10},{10,10}},
                  origin={-110,0})));

          Connectors.RealOutputBus y "Output bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50}), iconTransformation(extent={{-10,-10},{10,10}},
                  origin={110,0})));

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
        equation
          // Graphite
          connect(graphite.negative, negative.graphite) annotation (Line(
              points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));

          connect(graphite.positive, positive.graphite) annotation (Line(
              points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              pattern=LinePattern.None,
              thickness=0.5,
              smooth=Smooth.None));

          connect(u.graphite, graphite.u) annotation (Line(
              points={{5.55112e-16,50},{6.10623e-16,50},{6.10623e-16,5}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(graphite.y, y.graphite) annotation (Line(
              points={{6.10623e-16,-5},{5.55112e-16,-5},{5.55112e-16,-50}},
              color={0,0,127},
              thickness=0.5,
              smooth=Smooth.None));

          annotation (defaultComponentName="face", Diagram(graphics));
        end FlowsGraphiteOnly;

        model Efforts
          "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors, with differences in efforts specified by default</html>"
          extends Flows(
            gas(
              H2(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate),
              H2O(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate),
              N2(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate),
              O2(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate)),
            graphite('C+'(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate), 'e-'(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate)),
            ionomer(
              'SO3-'(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate),
              'H+'(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate),
              H2O(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate)),
            liquid(H2O(
                redeclare replaceable function materialSpec =
                    Face.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Pair.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Pair.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Pair.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Pair.Thermal.heatRate)));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          annotation (defaultComponentName="face",Diagram(graphics));
        end Efforts;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
          extends Modelica.Icons.Package;

          model Gas "Condition for gas"

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2=false "Include H2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows H2 if inclH2 "H2 conditions" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>H<sub>2</sub> conditions</html>",
                enable=inclH2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclN2=false "Include H2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Nitrogen (N<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows N2 if inclN2 "N2 conditions" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>N<sub>2</sub> conditions</html>",
                enable=inclN2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclO2=false "Include O2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Oxygen (O<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows O2 if inclO2 "O2 conditions" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>O<sub>2</sub> conditions</html>",
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclC+'=false "Include C+" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows 'C+' if 'inclC+' "C+ conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>C<sup>+</sup> conditions</html>",
                enable='inclC+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));
            // **Update

            parameter Boolean 'incle-'=false "Include e-" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows 'e-' if 'incle-' "e- conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>e<sup>-</sup> conditions</html>",
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclSO3-'=false
              "Include C19HF37O5S- (abbreviated as SO3-)" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label=
                    "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>, abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows 'SO3-' if 'inclSO3-' "SO3- conditions"
              annotation (Dialog(
                group="Species",
                __Dymola_label=
                    "<html>SO<sub>3</sub><sup>-</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclSO3-'), Placement(transformation(extent={{-10,-10},
                      {10,10}})));

            // **Update

            parameter Boolean 'inclH+'=false "Include H+" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows 'H+' if 'inclH+' "H+ conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sup>+</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclH+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // C19HF37O5S-
            connect('SO3-'.negative, negative.'SO3-') annotation (Line(
                points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},

                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('SO3-'.positive, positive.'SO3-') annotation (Line(
                points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(u.'SO3-', 'SO3-'.u) annotation (Line(
                points={{5.55112e-16,50},{6.10623e-16,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('SO3-'.y, y.'SO3-') annotation (Line(
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Pair.FaceFlows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
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

        protected
          model EmptyPhase "Empty condition for a phase (no species)"
            extends FCSys.Icons.Conditions.PairShort;

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

            Connectors.RealOutputBus y "Output bus of measurements" annotation
              (Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,-50}),iconTransformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,-50})));
            annotation (Icon(graphics));

          end EmptyPhase;
        end Phases;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"

        extends Modelica.Icons.Package;

        model Flows
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows specified by default</html>"

          extends FCSys.Icons.Conditions.SingleShort;

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
        end Flows;

        model FlowsFluidOnly
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows specified by default and only the fluid phases included</html>"

          extends FCSys.Icons.Conditions.SingleShort;

          Phases.Gas gas "Gas" annotation (Dialog(group=
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
        end FlowsFluidOnly;

        model FlowsGraphiteOnly
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows specified by default and only the graphite phase</html>"

          extends FCSys.Icons.Conditions.SingleShort;

          Phases.Graphite graphite "Graphite" annotation (Dialog(group=
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

          annotation (defaultComponentName="face");
        end FlowsGraphiteOnly;

        model Efforts
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with efforts specified by default</html>"
          extends Flows(
            gas(
              H2(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K)),
              H2O(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K)),
              N2(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K)),
              O2(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K))),
            graphite(redeclare replaceable ThermalDiffusion.Temperature 'C+',
                'e-'(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K))),
            ionomer(
              redeclare replaceable ThermalDiffusion.Temperature 'SO3-',
              'H+'(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K)),
              H2O(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K))),
            liquid(H2O(
                redeclare replaceable function materialSpec =
                    Face.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function beforeSpec =
                    Face.Single.Translational.velocity,
                redeclare replaceable function thermalSpec =
                    Face.Single.Thermal.temperature,
                redeclare replaceable function materialMeas =
                    Face.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function beforeMeas =
                    Face.Single.Translational.force,
                redeclare replaceable function thermalMeas =
                    Face.Single.Thermal.heatRate,
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  materialSet(y=U.atm),
                redeclare replaceable Modelica.Blocks.Sources.RealExpression
                  thermalSet(y=298.15*U.K))));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          // The daltonSource and thermalSet blocks are redeclared as not replaceable
          // because y is set directly and cannot be undone at instantiation.

          annotation (defaultComponentName="face",Diagram(graphics));
        end Efforts;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
          extends Modelica.Icons.Package;

          model Gas "Condition for gas"

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2=false "Include H2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows H2 if inclH2 "Include H2" annotation (Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclN2=false "Include N2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Nitrogen (N<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows N2 if inclN2 "N2 conditions" annotation (Dialog(
                group="Species",
                __Dymola_label="<html>N<sub>2</sub>Conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclN2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

            parameter Boolean inclO2=false "Include O2" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Oxygen (O<sub>2</sub>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows O2 if inclO2 "O2 conditions" annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>O<sub>2</sub> conditions</html>",
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclC+'=false "Include C+" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            replaceable ThermalDiffusion.HeatRate 'C+' if 'inclC+'
              constrainedby
              FCSys.Conditions.ByConnector.ThermalDiffusion.Partial
              "C+ conditions" annotation (
              __Dymola_choicesFromPackage=true,
              Dialog(
                group="Species",
                __Dymola_label="<html>C<sup>+</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclC+'),
              Placement(transformation(extent={{-10,-10},{10,10}})));

            parameter Boolean 'incle-'=false "Include e-" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows 'e-' if 'incle-' "e- conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>e<sup>-</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='incle-'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

          equation
            // C+
            connect('C+'.thermal, face.'C+') annotation (Line(
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean 'inclSO3-'=false
              "Iclude C19HF37O5S- (abbreviated as SO3-)" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label=
                    "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>, abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            replaceable ThermalDiffusion.HeatRate 'SO3-' if 'inclSO3-'
              constrainedby
              FCSys.Conditions.ByConnector.ThermalDiffusion.Partial
              "SO3- conditions" annotation (
              __Dymola_choicesFromPackage=true,
              Dialog(
                group="Species",
                __Dymola_label=
                    "<html>SO<sub>3</sub><sup>-</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclSO3-'),
              Placement(transformation(extent={{-10,-10},{10,10}})));

            parameter Boolean 'inclH+'=false "Include H+" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows 'H+' if 'inclH+' "H+ conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sup>+</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclH+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // C19HF37O5S-
            connect('SO3-'.thermal, face.'SO3-') annotation (Line(
                points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(u.'SO3-', 'SO3-'.u) annotation (Line(
                points={{-100,5.55112e-16},{-100,0},{-11,0},{-11,6.10623e-16}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('SO3-'.y, y.'SO3-') annotation (Line(
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

            extends EmptyPhase;

            // Conditionally include species.
            parameter Boolean inclH2O=false "Include H2O" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Face.Single.Flows H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
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

        protected
          model EmptyPhase "Empty condition for a phase (no species)"
            extends FCSys.Icons.Conditions.SingleShort;

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

            Connectors.RealOutputBus y "Output bus of measurements" annotation
              (Placement(transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={100,0}),iconTransformation(
                  extent={{-10,-10},{10,10}},
                  rotation=0,
                  origin={110,0})));
            annotation (Icon(graphics));

          end EmptyPhase;
        end Phases;

      end Single;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
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
        model FaceFlows
          "<html>Conditions for a <a href=\\\"modelica://FCSys.Connectors.Face\\\">Face</a> connector, with flows specified by default</html>"
          import Modelica.Blocks.Sources;
          extends FCSys.Icons.Conditions.PairShort;

          // Specification
          // -------------
          // Material
          replaceable function materialSpec = Material.current constrainedby
            Material.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            choicesAllMatching=true,
            Dialog(tab="Specification", group="Material"));
          parameter Boolean internalMaterial=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Material"));
          replaceable Sources.RealExpression materialSet if internalMaterial
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Material",
              enable=internalMaterial),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-50,40})));
          //
          // 1st transverse
          replaceable function afterSpec = Translational.force constrainedby
            Translational.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="First transverse"),
            Placement(transformation(extent={{-24,4},{-4,24}})));
          parameter Boolean internalAfter=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="First transverse"));
          replaceable Sources.RealExpression afterSet if internalAfter
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="First transverse",
              enable=internalAfter),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-10,40})));

          //
          // 2nd transverse
          replaceable function beforeSpec = Translational.force constrainedby
            Translational.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="Second transverse"),
            Placement(transformation(extent={{4,-10},{24,10}})));

          parameter Boolean internalBefore=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Second transverse"));
          replaceable Sources.RealExpression beforeSet if internalBefore
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Second transverse",
              enable=internalBefore),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={30,40})));

          //
          // Thermal
          replaceable function thermalSpec = Thermal.heatRate constrainedby
            Thermal.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="Thermal"),
            Placement(transformation(extent={{4,-10},{24,10}})));

          parameter Boolean internalThermal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Thermal"));
          replaceable Sources.RealExpression thermalSet if internalThermal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Thermal",
              enable=internalThermal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={70,40})));

          // Measurement
          // -----------
          // Material
          replaceable function materialMeas = Material.pressure constrainedby
            Material.Partial "Material quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // 1st transverse
          replaceable function afterMeas = Translational.velocity
            constrainedby Translational.Partial "First transverse quantity"
            annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                  "Measurement"));

          // 2nd transverse
          replaceable function beforeMeas = Translational.velocity
            constrainedby Translational.Partial "Second transverse quantity"
            annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                  "Measurement"));

          // Thermal
          replaceable function thermalMeas = Thermal.temperature constrainedby
            Thermal.Partial "Thermal quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // Aliases
          Q.Pressure Deltap "Difference in pressure";
          Q.Velocity Deltaphi[Orient] "Difference in velocity";
          Q.Temperature DeltaT "Difference in temperature";

          Connectors.Face negative "Negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}})));
          Connectors.Face positive "Positive face"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));
          Connectors.RealInputBus u if not (internalMaterial and internalAfter
             and internalBefore and internalThermal) "Bus of specifications"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,110}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50})));
          Connectors.RealOutputBus y "Bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-110}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50})));

          // Inputs
        protected
          Connectors.RealInputInternal u_material if not internalMaterial
            "Material specification" annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-60,70}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-80,110})));
          Connectors.RealInputInternal u_after if not internalAfter
            "First transverse specification" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-20,70}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,110})));
          Connectors.RealInputInternal u_before if not internalBefore
            "Second transverse specification" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={20,70}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={40,110})));
          Connectors.RealInputInternal u_thermal if not internalThermal
            "Thermal specification" annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,70}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={80,110})));

          // Outputs
          Connectors.RealOutputInternal _u_material=materialSpec(
                      Deltap,
                      negative.Ndot,
                      Deltaphi,
                      negative.mPhidot,
                      DeltaT,
                      negative.Qdot)
            "Internal, working value of material specification" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-60,6})));
          Connectors.RealOutputInternal _u_after=afterSpec(
                      Deltap,
                      negative.Ndot,
                      Deltaphi,
                      negative.mPhidot,
                      DeltaT,
                      negative.Qdot,
                      orient=Orient.after)
            "Internal, working value of first transverse specification"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-20,6})));
          Connectors.RealOutputInternal _u_before=beforeSpec(
                      Deltap,
                      negative.Ndot,
                      Deltaphi,
                      negative.mPhidot,
                      DeltaT,
                      negative.Qdot,
                      orient=Orient.before)
            "Internal, working value of second transverse specification"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={20,6})));
          Connectors.RealOutputInternal _u_thermal=thermalSpec(
                      Deltap,
                      negative.Ndot,
                      Deltaphi,
                      negative.mPhidot,
                      DeltaT,
                      negative.Qdot)
            "Internal, working value of thermal specification" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,6})));

          Sources.RealExpression materialOut(y=materialMeas(
                        Deltap,
                        negative.Ndot,
                        Deltaphi,
                        negative.mPhidot,
                        DeltaT,
                        negative.Qdot)) "Generate the material output"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-60,-70})));
          Sources.RealExpression beforeOut(y=beforeMeas(
                        Deltap,
                        negative.Ndot,
                        Deltaphi,
                        negative.mPhidot,
                        DeltaT,
                        negative.Qdot,
                        orient=Orient.before))
            "Generate the 2nd transverse output" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={20,-70})));
          Sources.RealExpression afterOut(y=afterMeas(
                        Deltap,
                        negative.Ndot,
                        Deltaphi,
                        negative.mPhidot,
                        DeltaT,
                        negative.Qdot,
                        orient=Orient.after))
            "Generate the 1st transverse output" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-20,-70})));
          Sources.RealExpression thermalOut(y=thermalMeas(
                        Deltap,
                        negative.Ndot,
                        Deltaphi,
                        negative.mPhidot,
                        DeltaT,
                        negative.Qdot)) "Generate the thermal output"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,-70})));
        equation
          // Differences in efforts
          Deltap = positive.p - negative.p;
          Deltaphi = positive.phi - negative.phi;
          DeltaT = positive.T - negative.T;

          // Conservation (without storage)
          0 = positive.Ndot + negative.Ndot "Material";
          {0,0} = positive.mPhidot + negative.mPhidot "Translational momentum";
          DeltaT = positive.Qdot + negative.Qdot "Energy";

          // Material
          connect(u_material, _u_material) annotation (Line(
              points={{-60,70},{-60,6}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(materialSet.y, _u_material) annotation (Line(
              points={{-50,29},{-50,20},{-60,20},{-60,6}},
              color={0,0,127},
              smooth=Smooth.None));

          // First transverse
          connect(u_after, _u_after) annotation (Line(
              points={{-20,70},{-20,6}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(afterSet.y, _u_after) annotation (Line(
              points={{-10,29},{-10,20},{-20,20},{-20,6}},
              color={0,0,127},
              smooth=Smooth.None));

          // Second transverse
          connect(u_before, _u_before) annotation (Line(
              points={{20,70},{20,6}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(beforeSet.y, _u_before) annotation (Line(
              points={{30,29},{30,20},{20,20},{20,6}},
              color={0,0,127},
              smooth=Smooth.None));

          // Thermal
          connect(u_thermal, _u_thermal) annotation (Line(
              points={{60,70},{60,6}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(thermalSet.y, _u_thermal) annotation (Line(
              points={{70,29},{70,20},{60,20},{60,6}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(materialOut.y, y.material) annotation (Line(
              points={{-60,-81},{-60,-90},{0,-90},{0,-110}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(afterOut.y, y.after) annotation (Line(
              points={{-20,-81},{-20,-90},{0,-90},{0,-110}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(beforeOut.y, y.before) annotation (Line(
              points={{20,-81},{20,-90},{0,-90},{0,-110}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(u_material, u.material) annotation (Line(
              points={{-60,70},{-60,90},{0,90},{0,110}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u_after, u.after) annotation (Line(
              points={{-20,70},{-20,80},{-20,80},{-20,90},{0,90},{0,110},{0,110}},

              color={0,0,127},
              smooth=Smooth.None));

          connect(u_before, u.before) annotation (Line(
              points={{20,70},{20,90},{-20,90},{-20,90},{0,90},{0,110},{0,110}},

              color={0,0,127},
              smooth=Smooth.None));

          connect(u_thermal, u.thermal) annotation (Line(
              points={{60,70},{60,90},{-20,90},{-20,90},{0,90},{0,110},{0,110}},

              color={0,0,127},
              smooth=Smooth.None));

          connect(thermalOut.y, y.thermal) annotation (Line(
              points={{60,-81},{60,-90},{0,-90},{0,-110}},
              color={0,0,127},
              smooth=Smooth.None));

          annotation (
            defaultComponentName="face",
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics),
            Icon(graphics));
        end FaceFlows;

        model FaceEfforts
          "<html>Conditions for a pair of <a href=\\\"modelica://FCSys.Connectors.Face\\\">Face</a> connectors, with difference in efforts specified by default</html>"
          extends FaceFlows(
            redeclare replaceable function materialSpec = Material.pressure,
            redeclare replaceable function afterSpec = Translational.velocity,
            redeclare replaceable function beforeSpec = Translational.velocity,

            redeclare replaceable function thermalSpec = Thermal.temperature,
            redeclare replaceable function materialMeas = Material.current,
            redeclare replaceable function afterMeas = Translational.force,
            redeclare replaceable function beforeMeas = Translational.force,
            redeclare replaceable function thermalMeas = Thermal.heatRate);

          // See note in ReactionEfforts.
          annotation (defaultComponentName="face");

        end FaceEfforts;

        package Material "Material conditions"
          extends Modelica.Icons.Package;

          function pressure "Difference in pressure"
            extends Partial;

          algorithm
            x := Deltap;
            annotation (Inline=true);
          end pressure;

          function current "Current"
            extends Partial;

          algorithm
            x := Ndot;
            annotation (Inline=true);
          end current;

          partial function Partial
            "Template of a function to select a material quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.Pressure Deltap
              "<html>Difference in pressure (&Delta;&<i>p</i>)</html>";
            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity Deltaphi[Orient]
              "<html>Difference in velocity (&Delta;&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute DeltaT
              "<html>Difference in temperature (&Delta;<i>T</i>)</html>";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Material;

        package Translational "Translational conditions"
          extends Modelica.Icons.Package;

          function velocity "Difference in velocity"
            extends Partial;

          algorithm
            x := Deltaphi[orient];
          end velocity;

          function force "Non-equilibrium force"
            extends Partial;

          algorithm
            x := mPhidot[orient];
          end force;

          partial function Partial
            "Template of a function to select a translational quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.Pressure Deltap
              "<html>Difference in pressure (&Delta;&<i>p</i>)</html>";

            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity Deltaphi[Orient]
              "<html>Difference in velocity (&Delta;&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute DeltaT
              "<html>Difference in temperature (&Delta;<i>T</i>)</html>";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            input Orient orient
              "Orientation of translational momentum w.r.t. the face";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Translational;

        package Thermal "Thermal conditions"
          extends Modelica.Icons.Package;

          function temperature "Difference in temperature"
            extends Partial;

          algorithm
            x := DeltaT;
          end temperature;

          function heatRate "Rate of thermal conduction"
            extends Partial;

          algorithm
            x := Qdot;
          end heatRate;

          partial function Partial
            "Template of a function to select a thermal quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.Pressure Deltap
              "<html>Difference in pressure (&Delta;&<i>p</i>)</html>";
            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity Deltaphi[Orient]
              "<html>Difference in velocity (&Delta;&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute DeltaT
              "<html>Difference in temperature (&Delta;<i>T</i>)</html>";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Thermal;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"
        extends Modelica.Icons.Package;
        model Flows
          "<html>Conditions for a <a href=\\\"modelica://FCSys.Connectors.Face\\\">Face</a> connector, with flows specified by default</html>"
          import Modelica.Blocks.Sources;
          extends FCSys.Icons.Conditions.SingleShort;

          // Specification
          // -------------
          // Material
          replaceable function materialSpec = Material.current constrainedby
            Material.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            choicesAllMatching=true,
            Dialog(tab="Specification", group="Material"));
          parameter Boolean internalMaterial=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Material"));
          replaceable Sources.RealExpression materialSet if internalMaterial
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Material",
              enable=internalMaterial),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-40,70})));
          //
          // 1st transverse
          replaceable function afterSpec = Translational.force constrainedby
            Translational.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="First transverse"),
            Placement(transformation(extent={{-24,4},{-4,24}})));
          parameter Boolean internalAfter=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="First transverse"));
          replaceable Sources.RealExpression afterSet if internalAfter
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="First transverse",
              enable=internalAfter),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-40,30})));

          //
          // 2nd transverse
          replaceable function beforeSpec = Translational.force constrainedby
            Translational.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="Second transverse"),
            Placement(transformation(extent={{4,-10},{24,10}})));

          parameter Boolean internalBefore=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Second transverse"));
          replaceable Sources.RealExpression beforeSet if internalBefore
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Second transverse",
              enable=internalBefore),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-40,-10})));

          //
          // Thermal
          replaceable function thermalSpec = Thermal.heatRate constrainedby
            Thermal.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="Thermal"),
            Placement(transformation(extent={{4,-10},{24,10}})));

          parameter Boolean internalThermal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Thermal"));
          replaceable Sources.RealExpression thermalSet if internalThermal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(
              tab="Specification",
              group="Thermal",
              enable=internalThermal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-40,-50})));

          // Measurement
          // -----------
          // Material
          replaceable function materialMeas = Material.density constrainedby
            Material.Partial "Material quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // 1st transverse
          replaceable function afterMeas = Translational.velocity
            constrainedby Translational.Partial "First transverse quantity"
            annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                  "Measurement"));

          // 2nd transverse
          replaceable function beforeMeas = Translational.velocity
            constrainedby Translational.Partial "Second transverse quantity"
            annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                  "Measurement"));

          // Thermal
          replaceable function thermalMeas = Thermal.temperature constrainedby
            Thermal.Partial "Thermal quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          Connectors.Face face
            "Connector to transport material, translational momentum, and thermal energy"
            annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
                iconTransformation(extent={{-10,-50},{10,-30}})));
          Connectors.RealInputBus u if not (internalMaterial and internalAfter
             and internalBefore and internalThermal) "Bus of specifications"
            annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
          Connectors.RealOutputBus y "Bus of measurements"
            annotation (Placement(transformation(extent={{100,-10},{120,10}})));

          // Inputs
        protected
          Connectors.RealInputInternal u_material if not internalMaterial
            "Material specification" annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,60})));
          Connectors.RealInputInternal u_after if not internalAfter
            "First transverse specification" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,20})));
          Connectors.RealInputInternal u_before if not internalBefore
            "Second transverse specification" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,-20})));
          Connectors.RealInputInternal u_thermal if not internalThermal
            "Thermal specification" annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-70,-60})));

          Connectors.RealOutputInternal _u_material
            "Internal, working value of material specification" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={4,60})));
          Connectors.RealOutputInternal _u_after=afterSpec(
                      face.p,
                      face.Ndot,
                      face.phi,
                      face.mPhidot,
                      face.T,
                      face.Qdot,
                      orient=Orient.after)
            "Internal, working value of first transverse specification"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={4,20})));
          Connectors.RealOutputInternal _u_before=beforeSpec(
                      face.p,
                      face.Ndot,
                      face.phi,
                      face.mPhidot,
                      face.T,
                      face.Qdot,
                      orient=Orient.before)
            "Internal, working value of second transverse specification"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={4,-20}), iconTransformation(extent={{-10,-10},{10,10}},
                  origin={6,-48})));
          Connectors.RealOutputInternal _u_thermal
            "Internal, working value of thermal specification" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={4,-60})));
        public
          Sources.RealExpression materialOut(y=materialMeas(
                        face.p,
                        face.Ndot,
                        face.phi,
                        face.mPhidot,
                        face.T,
                        face.Qdot)) "Generate the material output"
            annotation (Placement(transformation(extent={{40,50},{60,70}})));
          Sources.RealExpression beforeOut(y=beforeMeas(
                        face.p,
                        face.Ndot,
                        face.phi,
                        face.mPhidot,
                        face.T,
                        face.Qdot,
                        orient=Orient.before))
            "Generate the 2nd transverse output"
            annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
          Sources.RealExpression afterOut(y=afterMeas(
                        face.p,
                        face.Ndot,
                        face.phi,
                        face.mPhidot,
                        face.T,
                        face.Qdot,
                        orient=Orient.after))
            "Generate the 1st transverse output"
            annotation (Placement(transformation(extent={{40,10},{60,30}})));
          Sources.RealExpression thermalOut(y=thermalMeas(
                        face.p,
                        face.Ndot,
                        face.phi,
                        face.mPhidot,
                        face.T,
                        face.Qdot)) "Generate the thermal output"
            annotation (Placement(transformation(extent={{40,-70},{60,-50}})));
        equation
          _u_material = materialSpec(
                    face.p,
                    face.Ndot,
                    face.phi,
                    face.mPhidot,
                    face.T,
                    face.Qdot);
          _u_thermal = thermalSpec(
                    face.p,
                    face.Ndot,
                    face.phi,
                    face.mPhidot,
                    face.T,
                    face.Qdot);

          // Material
          connect(u_material, _u_material) annotation (Line(
              points={{-70,60},{4,60}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(materialSet.y, _u_material) annotation (Line(
              points={{-29,70},{-20,70},{-20,60},{4,60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u_material, u.material) annotation (Line(
              points={{-70,60},{-90,60},{-90,0},{-110,0}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(materialOut.y, y.material) annotation (Line(
              points={{61,60},{80,60},{80,0},{110,0}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{6,3},{6,3}}));

          // First transverse
          connect(u_after, _u_after) annotation (Line(
              points={{-70,20},{4,20}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(afterSet.y, _u_after) annotation (Line(
              points={{-29,30},{-20,30},{-20,20},{4,20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u_after, u.after) annotation (Line(
              points={{-70,20},{-90,20},{-90,0},{-110,0}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(afterOut.y, y.after) annotation (Line(
              points={{61,20},{80,20},{80,0},{84,0},{110,0}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{6,3},{6,3}}));

          // Second transverse
          connect(u_before, _u_before) annotation (Line(
              points={{-70,-20},{4,-20}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(beforeSet.y, _u_before) annotation (Line(
              points={{-29,-10},{-20,-10},{-20,-20},{4,-20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u_before, u.before) annotation (Line(
              points={{-70,-20},{-90,-20},{-90,0},{-110,0}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(beforeOut.y, y.before) annotation (Line(
              points={{61,-20},{80,-20},{80,0},{110,0}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{6,3},{6,3}}));

          // Thermal
          connect(thermalSet.y, _u_thermal) annotation (Line(
              points={{-29,-50},{-20,-50},{-20,-60},{4,-60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(u_thermal, u.thermal) annotation (Line(
              points={{-70,-60},{-90,-60},{-90,0},{-110,0}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(u_thermal, _u_thermal) annotation (Line(
              points={{-70,-60},{4,-60}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(thermalOut.y, y.thermal) annotation (Line(
              points={{61,-60},{80,-60},{80,0},{110,0}},
              color={0,0,127},
              smooth=Smooth.None), Text(
              string="%second",
              index=1,
              extent={{6,3},{6,3}}));
          annotation (
            defaultComponentName="face",
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics),
            Icon(graphics));
        end Flows;

        model Efforts
          "<html>Conditions for a pair of <a href=\\\"modelica://FCSys.Connectors.Face\\\">Face</a> connectors, with efforts specified by default</html>"
          extends Flows(
            redeclare replaceable function materialSpec = Material.pressure,
            redeclare replaceable function afterSpec = Translational.velocity,
            redeclare replaceable function beforeSpec = Translational.velocity,

            redeclare replaceable function thermalSpec = Thermal.temperature,
            redeclare replaceable function materialMeas = Material.current,
            redeclare replaceable function afterMeas = Translational.force,
            redeclare replaceable function beforeMeas = Translational.force,
            redeclare replaceable function thermalMeas = Thermal.heatRate,
            redeclare replaceable Modelica.Blocks.Sources.RealExpression
              materialSet(y=U.atm),
            redeclare replaceable Modelica.Blocks.Sources.RealExpression
              thermalSet(y=298.15*U.K));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          // The daltonSource and thermalSet blocks are redeclared as not replaceable
          // because y is set directly and cannot be undone at instantiation.

          annotation (defaultComponentName="face", Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));

        end Efforts;

        package Material "Material conditions"
          extends Modelica.Icons.Package;

          function pressure "Pressure"
            extends Partial;

          algorithm
            x := p;
            annotation (Inline=true);
          end pressure;

          function density "Density"
            extends Partial;

            replaceable package Data =
                Characteristics.BaseClasses.CharacteristicEOS constrainedby
              Characteristics.BaseClasses.CharacteristicEOS
              "Characteristic data" annotation (
              Dialog(group="Material properties"),
              choicesAllMatching=true,
              __Dymola_choicesFromPackage=true,
              Placement(transformation(extent={{-60,40},{-40,60}}),
                  iconTransformation(extent={{-10,90},{10,110}})));

          algorithm
            x := 1/Data.v_Tp(T, p);
            annotation (Inline=true);
          end density;

          function current "Current"
            extends Partial;

          algorithm
            x := Ndot;

            annotation (Inline=true);
          end current;

          function volumeRate "Volumic flow rate"
            extends Partial;

            replaceable package Data =
                Characteristics.BaseClasses.CharacteristicEOS constrainedby
              Characteristics.BaseClasses.CharacteristicEOS
              "Characteristic data" annotation (
              Dialog(group="Material properties"),
              choicesAllMatching=true,
              __Dymola_choicesFromPackage=true,
              Placement(transformation(extent={{-60,40},{-40,60}}),
                  iconTransformation(extent={{-10,90},{10,110}})));

          algorithm
            x := Data.v_Tp(T, p)*Ndot;
            annotation (Inline=true);
          end volumeRate;

          partial function Partial
            "Template of a function to select a material quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.PressureAbsolute p "<html>Pressure (<i>p</i>)</html>";
            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity phi[Orient] "<html>Velocity (&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute T "Temperature";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Material;

        package Translational "Translational conditions"
          extends Modelica.Icons.Package;

          function velocity "Velocity"
            extends Partial;

          algorithm
            x := phi[orient];
            annotation (Inline=true);
          end velocity;

          function force "Non-equilibrium force"
            extends Partial;

          algorithm
            x := mPhidot[orient];
            annotation (Inline=true);
          end force;

          partial function Partial
            "Template of a function to select a translational quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.PressureAbsolute p "<html>Pressure (<i>p</i>)</html>";
            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity phi[Orient] "<html>Velocity (&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute T "Temperature";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            input Orient orient
              "Orientation of translational momentum w.r.t. the face";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Translational;

        package Thermal "Thermal conditions"
          extends Modelica.Icons.Package;

          function temperature "Temperature"
            extends Partial;

          algorithm
            x := T;
            annotation (Inline=true);
          end temperature;

          function heatRate "Rate of thermal conduction"
            extends Partial;

          algorithm
            x := Qdot;
            annotation (Inline=true);
          end heatRate;

          partial function Partial
            "Template of a function to select a thermal quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.PressureAbsolute p "<html>Pressure (<i>p</i>)</html>";
            input Q.Current Ndot
              "<html>Diffusion current (<i>N&#775;</i>)</html>";

            // Translational
            input Q.Velocity phi[Orient] "<html>Velocity (&phi;)</html>";
            input Q.Force mPhidot[Orient]
              "<html>Non-equilibrium force (<i>m</i>&Phi;dot)</html>";

            // Thermal
            input Q.TemperatureAbsolute T "Temperature";
            input Q.Power Qdot
              "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

            output Real x "Value of condition";
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
          end Partial;
        end Thermal;

      end Single;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={191,191,191})}));

    end Face;

    package Amagat
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector</html>"
      extends Modelica.Icons.Package;

      model Pressure "Specify pressure (measure volume)"
        extends FCSys.Conditions.ByConnector.Amagat.Partial(final y=amagat.V,
            source(y=U.atm));

      equation
        amagat.p = u_final;

      end Pressure;

      model Volume "Specify volume (measure pressure)"
        extends Partial(final y=amagat.p, source(y=U.cc));

      equation
        amagat.V = u_final;

      end Volume;

      model Volume2 "Model to establish a fixed total volume"
        extends FCSys.Icons.Names.Top3;

        parameter Boolean setVolume=true
          "Set volume (otherwise, set pressure to the environment)"
          annotation (Dialog(compact=true), choices(__Dymola_checkBox=true));
        parameter Q.Volume V "Volume" annotation (Dialog(__Dymola_label=
                "<html><i>V</i></html>", enable=compressible));

        Connectors.Amagat amagat "Connector for additivity of volume"
          annotation (Placement(transformation(extent={{90,-50},{110,-30}}),
              iconTransformation(extent={{-10,-10},{10,10}})));

        outer Conditions.Environment environment "Environmental conditions";

      equation
        if setVolume then
          amagat.V = V;
        else
          amagat.p = environment.p;
        end if;

        annotation (
          defaultComponentName="volume",
          Documentation(info="<html><p>This model uses an <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector that imposes
    additivity of volume.  In order to use additivity of pressure, use
    the <a href=\"modelica://FCSys.Conditions.Adapters.AmagatDalton\">AmagatDalton</a> adapter.</p>

    <p>If <code>compressible</code> is <code>false</code>, the pressure is set (to the pressure
    of the <code>outer environment</code> model) instead of volume.</p>

    <p>See also the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={Polygon(
                      points={{-60,-60},{-60,20},{-20,60},{60,60},{60,-20},{20,
                  -60},{-60,-60}},
                      lineColor={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid)}),
          Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                  {100,100}}),graphics));
      end Volume2;

      partial model Partial "Base model for a pressure/volume"

        extends FCSys.Icons.Conditions.SingleShort;

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

        Connectors.RealOutput y "Measurement expression" annotation (Dialog(tab
              ="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.Amagat amagat "Connector for additivity of volume"
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
            points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{-20,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentName="amagat", Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics));
      end Partial;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}), Ellipse(
              extent={{-30,30},{30,-30}},
              fillColor={47,107,251},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0})}));
    end Amagat;

    package Dalton
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector</html>"
      extends Modelica.Icons.Package;

      model Volume "Specify volume (measure pressure)"
        extends Partial(final y=dalton.p, source(y=U.cc));

      equation
        dalton.V = u_final;

      end Volume;

      model Pressure "Specify pressure (measure volume)"
        extends Partial(final y=dalton.V, source(y=U.atm));

      equation
        dalton.p = u_final;

      end Pressure;

      partial model Partial "Base model for a pressure/volume"

        extends FCSys.Icons.Conditions.SingleShort;

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

        Connectors.RealOutput y "Measurement expression" annotation (Dialog(tab
              ="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        Connectors.Dalton dalton "Connector for additivity of pressure"
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
            points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{-20,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentName="dalton", Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics));
      end Partial;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={47,107,251})}));
    end Dalton;

    package Direct
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Direct\">Direct</a> or <a href=\"modelica://FCSys.Connectors.DirectNote\">DirectNote</a> connector</html>"
      extends Modelica.Icons.Package;

      model Flows
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Direct\">Direct</a> connector, with flows specified by default</html>"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.enumerate;
        import Modelica.Blocks.Sources;
        extends FCSys.Icons.Conditions.Single;

        // Specification
        // -------------
        // X-axis translational
        replaceable function transXSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX),
          Placement(transformation(extent={{-52,18},{-32,38}})));

        parameter Boolean internalTransX=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX));
        replaceable Sources.RealExpression transXSet if inclTransX and
          internalTransX constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX and internalTransX),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,70})));

        //
        // Y-axis translational
        replaceable function transYSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY),
          Placement(transformation(extent={{-24,4},{-4,24}})));

        parameter Boolean internalTransY=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY));
        replaceable Sources.RealExpression transYSet if inclTransY and
          internalTransY constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY and internalTransY),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,30})));
        //
        // Z-axis translational
        replaceable function transZSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalTransZ=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ));
        replaceable Sources.RealExpression transZSet if inclTransZ and
          internalTransZ constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ and internalTransZ),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-10})));
        //
        // Thermal
        replaceable function thermalSpec = Thermal.heatRate constrainedby
          Thermal.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Specification", group="Thermal"),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalThermal=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(tab="Specification", group="Thermal"));
        replaceable Sources.RealExpression thermalSet if internalThermal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Thermal",
            enable=internalThermal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-50})));

        // Measurement
        // -----------
        // X-axis translational
        replaceable function transXMeas = Translational.velocity constrainedby
          Translational.Partial "X-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Y-axis translational
        replaceable function transYMeas = Translational.velocity constrainedby
          Translational.Partial "Y-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Z-axis translational
        replaceable function transZMeas = Translational.velocity constrainedby
          Translational.Partial "Z-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Thermal
        replaceable function thermalMeas = Thermal.temperature constrainedby
          Thermal.Partial "Thermal quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));

        FCSys.Connectors.Direct direct(final n_trans=n_trans)
          "Connector to exchange translational momentum and thermal energy by diffusion"
          annotation (choicesAllMatching=true, Placement(transformation(extent=
                  {{-10,-110},{10,-90}})));

        // Inputs
        Connectors.RealInput u_transX if inclTransX and not internalTransX
          "X-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,60})));
        Connectors.RealInput u_transY if inclTransY and not internalTransY
          "Y-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,20})));
        Connectors.RealInput u_transZ if inclTransZ and not internalTransZ
          "Z-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-20})));
        Connectors.RealInput u_thermal if not internalThermal
          "Thermal specification" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-60})));

        // Outputs
        final Connectors.RealOutput y_transX=transXMeas(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,60}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,20})));

        final Connectors.RealOutput y_transY=transYMeas(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,20}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,60})));

        final Connectors.RealOutput y_transZ=transZMeas(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-20}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-20})));

        final Connectors.RealOutput y_thermal=thermalMeas(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot) "Thermal measurement" annotation (Dialog(
              tab="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-60})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

        Connectors.RealOutputInternal _u_transX=transXSpec(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "Internal, working value of X-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,60})));

        Connectors.RealOutputInternal _u_transY=transYSpec(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Internal, working value of Y-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,20})));

        Connectors.RealOutputInternal _u_transZ=transZSpec(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Internal, working value of Z-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-20})));

        Connectors.RealOutputInternal _u_thermal=thermalSpec(
                  direct.translational.phi,
                  direct.translational.mPhidot,
                  direct.thermal.T,
                  direct.thermal.Qdot)
          "Internal, working value of thermal specification" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-60})));
      equation
        // X-axis translational
        connect(u_transX, _u_transX) annotation (Line(
            points={{-110,60},{-36,60}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transXSet.y, _u_transX) annotation (Line(
            points={{-69,70},{-60,70},{-60,60},{-36,60}},
            color={0,0,127},
            smooth=Smooth.None));

        // Y-axis translational
        connect(u_transY, _u_transY) annotation (Line(
            points={{-110,20},{-36,20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transYSet.y, _u_transY) annotation (Line(
            points={{-69,30},{-60,30},{-60,20},{-36,20}},
            color={0,0,127},
            smooth=Smooth.None));

        // Z-axis translational
        connect(u_transZ, _u_transZ) annotation (Line(
            points={{-110,-20},{-36,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transZSet.y, _u_transZ) annotation (Line(
            points={{-69,-10},{-60,-10},{-60,-20},{-36,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        // Thermal
        connect(u_thermal, _u_thermal) annotation (Line(
            points={{-110,-60},{-36,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(thermalSet.y, _u_thermal) annotation (Line(
            points={{-69,-50},{-60,-50},{-60,-60},{-36,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (
          defaultComponentName="direct",
          Diagram(graphics),
          Icon(graphics));
      end Flows;

      model Efforts
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Direct\">Direct</a> connector, with efforts specified by default</html>"

        extends Flows(
          redeclare replaceable function transXSpec = Translational.velocity,
          redeclare replaceable function transYSpec = Translational.velocity,
          redeclare replaceable function transZSpec = Translational.velocity,
          redeclare replaceable function thermalSpec = Thermal.temperature,
          redeclare replaceable function transXMeas = Translational.force,
          redeclare replaceable function transYMeas = Translational.force,
          redeclare replaceable function transZMeas = Translational.force,
          redeclare replaceable function thermalMeas = Thermal.heatRate,
          redeclare replaceable Modelica.Blocks.Sources.RealExpression
            thermalSet(y=298.15*U.K));

        // The daltonSource and thermalSet blocks are redeclared as not replaceable
        // because y is set directly and cannot be undone at instantiation.

        annotation (defaultComponentName="direct");

      end Efforts;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;

        function velocity "Velocity"
          extends Partial;

        algorithm
          x := phi[i];
          annotation (Inline=true);
        end velocity;

        function force "Force"
          extends Partial;

        algorithm
          x := mPhidot[i];
          annotation (Inline=true);
        end force;

        partial function Partial
          "Template of a function to select a translational quantity"
          extends Modelica.Icons.Function;

          // Translational
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal
          input Q.TemperatureAbsolute T "Temperature";
          input Q.Power Qdot
            "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

          input Integer i(min=1,max=3) "Index of the translational axis";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Translational;

      package Thermal "Thermal conditions"
        extends Modelica.Icons.Package;

        function temperature "Temperature"
          extends Partial;

        algorithm
          x := T;
          annotation (Inline=true);
        end temperature;

        function heatRate "Heat flow rate"
          extends Partial;

        algorithm
          x := Qdot;
          annotation (Inline=true);
        end heatRate;

        partial function Partial
          "Template of a function to select a thermal quantity"
          extends Modelica.Icons.Function;

          // Translational
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal
          input Q.TemperatureAbsolute T "Temperature";
          input Q.Power Qdot
            "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Thermal;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={2,157,21},
              fillPattern=FillPattern.Solid,
              fillColor={38,196,52})}));
    end Direct;

    package Inert
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>, <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a>, or <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> connector</html>"
      model Flows
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a> or <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> connector, with flows specified by default</html>"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.enumerate;
        import Modelica.Blocks.Sources;
        extends FCSys.Icons.Conditions.Single;

        // Specification
        // -------------
        // X-axis translational
        replaceable function transXSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX),
          Placement(transformation(extent={{-52,18},{-32,38}})));

        parameter Boolean internalTransX=true if inclTransX
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX));
        replaceable Sources.RealExpression transXSet if inclTransX and
          internalTransX constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX and internalTransX),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,70})));

        //
        // Y-axis translational
        replaceable function transYSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY),
          Placement(transformation(extent={{-24,4},{-4,24}})));

        parameter Boolean internalTransY=true if inclTransY
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY));
        replaceable Sources.RealExpression transYSet if inclTransY and
          internalTransY constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY and internalTransY),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,30})));
        //
        // Z-axis translational
        replaceable function transZSpec = Translational.force constrainedby
          Translational.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalTransZ=true if inclTransZ
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ));
        replaceable Sources.RealExpression transZSet if inclTransZ and
          internalTransZ constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ and internalTransZ),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-10})));
        //
        // Thermal
        replaceable function thermalSpec = Thermal.heatRate constrainedby
          Thermal.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Specification", group="Thermal"),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalThermal=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(tab="Specification", group="Thermal"));
        replaceable Sources.RealExpression thermalSet if internalThermal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Thermal",
            enable=internalThermal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-50})));

        // Measurement
        // -----------
        // X-axis translational
        replaceable function transXMeas =
            FCSys.Conditions.ByConnector.Inert.Translational.velocity
          constrainedby Translational.Partial "X-axis translational quantity"
          annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                "Measurement"));
        //
        // Y-axis translational
        replaceable function transYMeas =
            FCSys.Conditions.ByConnector.Inert.Translational.velocity
          constrainedby Translational.Partial "Y-axis translational quantity"
          annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                "Measurement"));
        //
        // Z-axis translational
        replaceable function transZMeas =
            FCSys.Conditions.ByConnector.Inert.Translational.velocity
          constrainedby Translational.Partial "Z-axis translational quantity"
          annotation (__Dymola_choicesFromPackage=true, Dialog(tab=
                "Measurement"));
        //
        // Thermal
        replaceable function thermalMeas =
            FCSys.Conditions.ByConnector.Inert.Thermal.temperature
          constrainedby Thermal.Partial "Thermal quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));

        Connectors.Inter inter(final n_trans=n_trans)
          "Connector to exchange translational momentum and thermal energy by diffusion"
          annotation (choicesAllMatching=true, Placement(transformation(extent=
                  {{-10,-110},{10,-90}})));

        // Inputs
        Connectors.RealInput u_transX if inclTransX and not internalTransX
          "X-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,60})));
        Connectors.RealInput u_transY if inclTransY and not internalTransY
          "Y-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,20})));
        Connectors.RealInput u_transZ if inclTransZ and not internalTransZ
          "Z-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-20})));
        Connectors.RealInput u_thermal if not internalThermal
          "Thermal specification" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-60})));

        // Outputs
        final Connectors.RealOutput y_transX=transXMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,60}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        final Connectors.RealOutput y_transY=transYMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,20}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));

        final Connectors.RealOutput y_transZ=transZMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-20}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40})));

        final Connectors.RealOutput y_thermal=thermalMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot) "Thermal measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-80})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

        Connectors.RealOutputInternal _u_transX=transXSpec(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "Internal, working value of X-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,60})));

        Connectors.RealOutputInternal _u_transY=transYSpec(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Internal, working value of Y-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,20})));

        Connectors.RealOutputInternal _u_transZ=transZSpec(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Internal, working value of Z-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-20})));

        Connectors.RealOutputInternal _u_thermal=thermalSpec(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot)
          "Internal, working value of thermal specification" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-60})));

      equation
        // X-axis translational
        connect(u_transX, _u_transX) annotation (Line(
            points={{-110,60},{-36,60}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transXSet.y, _u_transX) annotation (Line(
            points={{-69,70},{-60,70},{-60,60},{-36,60}},
            color={0,0,127},
            smooth=Smooth.None));

        // Y-axis translational
        connect(u_transY, _u_transY) annotation (Line(
            points={{-110,20},{-36,20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transYSet.y, _u_transY) annotation (Line(
            points={{-69,30},{-60,30},{-60,20},{-36,20}},
            color={0,0,127},
            smooth=Smooth.None));

        // Z-axis translational
        connect(u_transZ, _u_transZ) annotation (Line(
            points={{-110,-20},{-36,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transZSet.y, _u_transZ) annotation (Line(
            points={{-69,-10},{-60,-10},{-60,-20},{-36,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        // Thermal
        connect(u_thermal, _u_thermal) annotation (Line(
            points={{-110,-60},{-36,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(thermalSet.y, _u_thermal) annotation (Line(
            points={{-69,-50},{-60,-50},{-60,-60},{-36,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (
          defaultComponentName="inert",
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics),
          Icon(graphics));
      end Flows;
      extends Modelica.Icons.Package;
      model Efforts
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a> or <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> connector, with efforts specified by default</html>"

        extends FCSys.Conditions.ByConnector.Inert.Flows(
          redeclare replaceable function transXSpec = Translational.velocity,
          redeclare replaceable function transYSpec = Translational.velocity,
          redeclare replaceable function transZSpec = Translational.velocity,
          redeclare replaceable function thermalSpec = Thermal.temperature,
          redeclare replaceable function transXMeas = Translational.force,
          redeclare replaceable function transYMeas = Translational.force,
          redeclare replaceable function transZMeas = Translational.force,
          redeclare replaceable function thermalMeas = Thermal.heatRate,
          redeclare replaceable Modelica.Blocks.Sources.RealExpression
            thermalSet(y=298.15*U.K));

        // The thermalSet block is redeclared as not replaceable
        // because y is set directly and cannot be undone at instantiation.

        annotation (defaultComponentName="inert");

      end Efforts;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;

        function velocity "Velocity"
          extends Partial;

        algorithm
          x := phi[i];
          annotation (Inline=true);
        end velocity;

        function force "Force"
          extends Partial;

        algorithm
          x := mPhidot[i];
          annotation (Inline=true);
        end force;

        partial function Partial
          "Template of a function to select a translational quantity"
          extends Modelica.Icons.Function;

          // Translational
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal
          input Q.TemperatureAbsolute T "Temperature";
          input Q.Power Qdot
            "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

          input Integer i(min=1,max=3) "Index of the translational axis";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Translational;

      package Thermal "Thermal conditions"
        extends Modelica.Icons.Package;

        function temperature "Temperature"
          extends Partial;

        algorithm
          x := T;
          annotation (Inline=true);
        end temperature;

        function heatRate "Heat flow rate"
          extends Partial;

        algorithm
          x := Qdot;
          annotation (Inline=true);
        end heatRate;

        partial function Partial
          "Template of a function to select a thermal quantity"
          extends Modelica.Icons.Function;

          // Translational
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          // Thermal
          input Q.TemperatureAbsolute T "Temperature";
          input Q.Power Qdot
            "<html>Rate of thermal conduction (<i>Q&#775;</i>)</html>";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Thermal;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={2,157,21},
              fillPattern=FillPattern.Solid,
              fillColor={38,196,52})}));
    end Inert;

    package Translational
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector</html>"
      extends Modelica.Icons.Package;

      model Force
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector, with force specified by default</html>"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.enumerate;
        import Modelica.Blocks.Sources;
        extends FCSys.Icons.Conditions.Single;

        // Specification
        // -------------
        // X-axis translational
        replaceable function transXSpec = Component.force constrainedby
          Component.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX),
          Placement(transformation(extent={{-52,18},{-32,38}})));

        parameter Boolean internalTransX=true if inclTransX
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX));
        replaceable Sources.RealExpression transXSet if inclTransX and
          internalTransX constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="X-axis translational",
            enable=inclTransX and internalTransX),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,50})));

        //
        // Y-axis translational
        replaceable function transYSpec = Component.force constrainedby
          Component.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY),
          Placement(transformation(extent={{-24,4},{-4,24}})));

        parameter Boolean internalTransY=true if inclTransY
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY));
        replaceable Sources.RealExpression transYSet if inclTransY and
          internalTransY constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Y-axis translational",
            enable=inclTransY and internalTransY),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,10})));

        //
        // Z-axis translational
        replaceable function transZSpec = Component.force constrainedby
          Component.Partial "Quantity" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ),
          Placement(transformation(extent={{4,-10},{24,10}})));

        parameter Boolean internalTransZ=true if inclTransZ
          "Use internal specification" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ));
        replaceable Sources.RealExpression transZSet if inclTransZ and
          internalTransZ constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            tab="Specification",
            group="Z-axis translational",
            enable=inclTransZ and internalTransZ),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-80,-30})));

        // Measurement
        // -----------
        // X-axis translational
        replaceable function transXMeas = Component.velocity constrainedby
          Component.Partial "X-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Y-axis translational
        replaceable function transYMeas = Component.velocity constrainedby
          Component.Partial "Y-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Z-axis translational
        replaceable function transZMeas = Component.velocity constrainedby
          Component.Partial "Z-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Included components of translational momentum
        parameter Boolean inclTransX=true "X" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransY=true "Y" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));
        parameter Boolean inclTransZ=true "Z" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            tab="Assumptions",
            group="Included transport axes",
            compact=true));

        // Inputs
        Connectors.RealInput u_transX if inclTransX and not internalTransX
          "X-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,40})));
        Connectors.RealInput u_transY if inclTransY and not internalTransY
          "Y-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,0})));
        Connectors.RealInput u_transZ if inclTransZ and not internalTransZ
          "Z-axis translational specification" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-110,-40})));

        // Outputs
        final Connectors.RealOutput y_transX=transXMeas(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        final Connectors.RealOutput y_transY=transYMeas(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));

        final Connectors.RealOutput y_transZ=transZMeas(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Dialog(group=
                "Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40})));

        Connectors.Translational translational(final n_trans=n_trans)
          "Connector for advection or diffusion of translational momentum"
          annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
              iconTransformation(extent={{-10,-110},{10,-90}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

        Connectors.RealOutputInternal _u_transX=transXSpec(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.x]) if inclTransX
          "Internal, working value of X-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,40})));

        Connectors.RealOutputInternal _u_transY=transYSpec(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.y]) if inclTransY
          "Internal, working value of Y-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,0})));

        Connectors.RealOutputInternal _u_transZ=transZSpec(
                  translational.phi,
                  translational.mPhidot,
                  i=transCart[Axis.z]) if inclTransZ
          "Internal, working value of Z-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,-40})));

      equation
        // X-axis translational
        connect(u_transX, _u_transX) annotation (Line(
            points={{-110,40},{-36,40}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transXSet.y, _u_transX) annotation (Line(
            points={{-69,50},{-60,50},{-60,40},{-36,40}},
            color={0,0,127},
            smooth=Smooth.None));

        // Y-axis translational
        connect(u_transY, _u_transY) annotation (Line(
            points={{-110,5.55112e-16},{-36,5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transYSet.y, _u_transY) annotation (Line(
            points={{-69,10},{-60,10},{-60,5.55112e-16},{-36,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        // Z-axis translational
        connect(u_transZ, _u_transZ) annotation (Line(
            points={{-110,-40},{-36,-40}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(transZSet.y, _u_transZ) annotation (Line(
            points={{-69,-30},{-60,-30},{-60,-40},{-36,-40}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (Diagram(graphics), Icon(graphics));
      end Force;

      model Velocity
        "<html>Condition for a <a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector, with velocity specified by default</html>"

        extends Force(
          redeclare replaceable Component.velocity transXSpec,
          redeclare replaceable Component.velocity transYSpec,
          redeclare replaceable Component.velocity transZSpec,
          redeclare replaceable Component.force transXMeas,
          redeclare replaceable Component.force transYMeas,
          redeclare replaceable Component.force transZMeas);

      end Velocity;

      package Component "Conditions for a component of translational momentum"
        extends Modelica.Icons.Package;
        function velocity "Velocity "
          extends Partial;

        algorithm
          x := phi[i];
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="translational");
        end velocity;

        function force "Force"
          extends Partial;

        algorithm
          x := mPhidot[i];

        end force;

        partial function Partial
          "Template of a function to select a translational quantity"
          extends Modelica.Icons.Function;

          // Translational advection
          input Q.Velocity phi[:] "<html>Velocity (&phi;)</html>";
          input Q.Force mPhidot[:] "<html>Force (<i>m</i>&Phi;dot)</html>";

          input Integer i(min=1,max=3) "Index of the translational axis";

          output Real x "Value of condition";
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));
        end Partial;
      end Component;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255})}));

    end Translational;

    package ThermalDiffusion
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.ThermalDiffusion\">ThermalDiffusion</a> connector</html>"
      extends Modelica.Icons.Package;

      model Temperature "Specify temperature (measure heat flow rate)"
        extends Partial(final y=thermal.Qdot, source(y=298.15*U.K));

      equation
        thermal.T = u_final;

      end Temperature;

      model HeatRate "Specify heat flow rate (measure temperature)"
        extends Partial(final y=thermal.T);

      equation
        thermal.Qdot = u_final;

      end HeatRate;

      partial model Partial "Base model for a thermal condition"

        extends FCSys.Icons.Conditions.SingleShort;

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

        Connectors.RealOutput y "Measurement expression" annotation (Dialog(tab
              ="Measurement"), Placement(transformation(
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
            points={{-110,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{-20,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentName="thermal");
      end Partial;
      annotation (Icon(graphics={Ellipse(
              extent={{-60,60},{60,-60}},
              lineColor={127,127,127},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255})}));
    end ThermalDiffusion;
    annotation (Documentation(info="<html>
  <p>This package contains models to impose conditions on each of the declarative connectors
  established in <a href=\"modelica://FCSys.Connectors\">FCSys.Connectors</a>.  The subpackages
  are named according to the corresponding connector.</p>
</html>"));
  end ByConnector;

  package TestStands "Test stands"
    extends Modelica.Icons.Package;
    extends Modelica.Icons.UnderConstruction;

    model TestStandEIS
      "Test stand to perform electrochemical impedance spectroscopy"
      extends TestStand(redeclare Q.Current zI,zJ=zJ_large + zJ_small_SI*U.A/U.m
            ^2);

      parameter Q.CurrentAreic zJ_large=U.A "Large-signal current density"
        annotation (Dialog(__Dymola_label=
              "<html><i>zJ</i><sub>large</sub></html>"));
      Connectors.RealInput zJ_small_SI
        "Small-signal current density in SI base units" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-45,
            origin={-107,107}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=-45,
            origin={-167,167})));
      Connectors.RealOutput w_V "Cell potential in volts" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=-45,
            origin={107,-107}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=-45,
            origin={167,-167})));

    equation
      w_V = w/U.V;

      annotation (
        Documentation(info="<html><p>This model modulates the electrical current applied to the cell 
    according to an input.
    The current density is the sum of a steady-state large-signal current density and a small-signal 
    current density introduced via the input <i>zJ</i><sub>small SI</sub>.</p>
       
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Conditions.TestStands.TestStand\">test stand</a> model.</p></html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
                100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{
                160,160}}), graphics));
    end TestStandEIS;

    model TestStand "Fuel cell test stand (applies boundary conditions)"
      import saturationPressureSI =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid;
      import FCSys.Utilities.average;
      import FCSys.Utilities.inSign;
      import FCSys.Conditions.ByConnector.Face.Single;
      extends FCSys.Icons.Names.Top9;

      // Geometry
      parameter Q.Length L_x_an[:]={8*U.mm}
        "Lengths of the segments through the cell in anode FP" annotation (
          Dialog(group="Cell geometry", __Dymola_label=
              "<html><i>L</i><sub>x an</sub></html>"));
      parameter Q.Length L_x_ca[:]={8*U.mm}
        "Lengths of the segments through the cell in cathode FP" annotation (
          Dialog(group="Cell geometry", __Dymola_label=
              "<html><i>L</i><sub>x ca</sub></html>"));
      parameter Q.Length L_y[:]={U.m}
        "Lengths of the segments along the channel" annotation (Dialog(group=
              "Cell geometry", __Dymola_label=
              "<html><i>L</i><sub>y</sub></html>"));
      parameter Q.Length L_z[:]={5*U.mm}
        "Lengths of the segments across the channel" annotation (Dialog(group=
              "Cell geometry",__Dymola_label=
              "<html><i>L</i><sub>z</sub></html>"));
      final parameter Integer n_x_an=size(L_x_an, 1)
        "Number of subregions along the through-cell axis in anode FP"
        annotation (Dialog(group="Cell geometry"));
      final parameter Integer n_x_ca=size(L_x_ca, 1)
        "Number of subregions along the through-cell axis in cathode FP"
        annotation (Dialog(group="Cell geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";
      final parameter Q.Area A=sum(L_y)*sum(L_z) "Cross-sectional area";
      final parameter Q.Area A_seg[n_y, n_z]=outerProduct(L_y, L_z)
        "Areas of the yz segments";
      final parameter Q.Area A_an_seg[n_x_an, n_z]=outerProduct(L_x_an, L_z)
        "Areas of the xz segments of the anode flow plate";
      final parameter Q.Area A_ca_seg[n_x_ca, n_z]=outerProduct(L_x_ca, L_z)
        "Areas of the xz segments of the cathode flow plate";
      final parameter Q.Area A_an=sum(L_x_an)*sum(L_z)
        "Total cross-sectional area of the anode flow plate in the xz plane";
      final parameter Q.Area A_ca=sum(L_x_ca)*sum(L_z)
        "Total cross-sectional area of the cathode flow plate in the xz plane";

      // Operating conditions
      // --------------------
      // Electrical
      parameter Enumerations.ElectricalSpec electricalSpec=ElectricalSpec.currentDensity
        "Type of electrical specification" annotation (Dialog(
          tab="Conditions",
          group="Electrical",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of specification",
          __Dymola_joinNext=true));
      Real u_electrical=U.A/U.cm^2 "Value of the electrical specification"
        annotation (Dialog(
          tab="Conditions",
          group="Electrical",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>electrical</sub>)</html>"));
      Q.CurrentAreic zJ "Current density";
      Q.Current zI "Current";
      Q.Potential w "Voltage";
      Q.ResistanceElectrical R "Resistance";
      Q.Power P "Power";
      //
      // General anode conditions
      parameter Side anInletSide=Side.p "Side of the inlet"
        annotation (Dialog(tab="Conditions",group="Anode"));
      Q.TemperatureAbsolute T_an_in=333.15*U.K "Inlet temperature" annotation (
          Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_label="<html><i>T</i><sub>an in</sub></html>"));

      Q.PressureAbsolute p_an_out=U.from_kPag(48.3) "Outlet pressure"
        annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_label="<html><i>p</i><sub>an out</sub></html>"));
      //
      // General cathode conditions
      parameter Side caInletSide=Side.p "Side of the inlet"
        annotation (Dialog(tab="Conditions",group="Cathode"));
      Q.TemperatureAbsolute T_ca_in=333.15*U.K "Inlet temperature" annotation (
          Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_label="<html><i>T</i><sub>ca in</sub></html>"));
      Q.PressureAbsolute p_ca_out=U.from_kPag(48.3) "Outlet pressure"
        annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_label="<html><i>p</i><sub>ca out</sub></html>"));
      Q.NumberAbsolute n_O2_in(
        final max=1,
        displayUnit="%") = 0.208
        "<html>Dry-gas concentration of O<sub>2</sub> at the inlet</html>"
        annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_label="<html><i>n</i><sub>O2 in</sub></html>"));
      //
      // Anode flow rate
      parameter FCSys.Conditions.TestStands.Enumerations.FlowSpec anFlowSpec=
          FlowSpec.stoich "Type of anode flow specification" annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of flow specification",
          __Dymola_joinNext=true));
      Real u_an_flow=1.5 "Value of the anode flow specification" annotation (
          Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>an flow</sub>)</html>"));
      Q.NumberAbsolute anStoich "Anode stoichiometric flow rate";
      Q.CurrentAreic J_an "Equivalent current density of anode supply";
      Q.Current I_an "Equivalent current of anode supply";
      Q.VolumeRate Vdot_g_an_in "Volumetric flow rate of gas in anode supply";
      Q.PressureAbsolute p_an_in "Anode inlet pressure";
      //
      // Cathode flow rate
      parameter FCSys.Conditions.TestStands.Enumerations.FlowSpec caFlowSpec=
          FlowSpec.stoich "Type of cathode flow specification" annotation (
          Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of flow specification",
          __Dymola_joinNext=true));
      Real u_ca_flow=2.0 "Value of the cathode flow specification" annotation (
          Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>ca flow</sub>)</html>"));
      Q.NumberAbsolute caStoich "Cathode stoichiometric flow rate";
      Q.CurrentAreic J_ca "Equivalent current density of cathode supply";
      Q.Current I_ca "Equivalent current of cathode supply";
      Q.VolumeRate Vdot_g_ca_in "Volumetric flow rate of gas in cathode supply";
      Q.PressureAbsolute p_ca_in "Cathode inlet pressure";
      //
      // Anode humidity
      parameter FCSys.Conditions.TestStands.Enumerations.HumiditySpec
        anHumiditySpec=HumiditySpec.relative
        "Type of anode humidity specification" annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of humidity specification",
          __Dymola_joinNext=true));
      Real u_an_humidity=0.8 "Value of the anode humidity specification"
        annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>an humidity</sub>)</html>"));
      Q.NumberAbsolute anInletRH(displayUnit="%")
        "Relative humidity at anode inlet";
      Q.PressureAbsolute p_H2O_an_in "H2O vapor pressure at anode inlet";
      Q.TemperatureAbsolute T_sat_an_in "Dew point at anode inlet";
      //
      // Cathode humidity
      parameter FCSys.Conditions.TestStands.Enumerations.HumiditySpec
        caHumiditySpec=HumiditySpec.relative
        "Type of anode humidity specification" annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of humidity specification",
          __Dymola_joinNext=true));
      Real u_ca_humidity=0.5 "Value of the cathode humidity specification"
        annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>ca humidity</sub>)</html>"));
      Q.NumberAbsolute caInletRH(displayUnit="%")
        "Relative humidity at cathode inlet";
      Q.PressureAbsolute p_H2O_ca_in "H2O vapor pressure at cathode inlet";
      Q.TemperatureAbsolute T_sat_ca_in "Dew point at cathode inlet";
      //
      // Anode end plate
      parameter FCSys.Conditions.TestStands.Enumerations.ThermalSpec
        anEndPlateSpec=ThermalSpec.temperature
        "Type of anode end plate specification" annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of end plate specification",
          __Dymola_joinNext=true));
      Real u_an_end_plate=333.15*U.K
        "Value of the anode end plate specification" annotation (Dialog(
          tab="Conditions",
          group="Anode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>an end plate</sub>)</html>"));
      Q.TemperatureAbsolute T_an "Temperature of anode end plate";
      Q.Conductance G_an
        "Thermal conductance of the anode end plate to the environment";
      Q.Power Qdot_an "Heat flow rate from the anode end plate";
      //
      // Cathode end plate
      parameter FCSys.Conditions.TestStands.Enumerations.ThermalSpec
        caEndPlateSpec=ThermalSpec.temperature
        "Type of anode end plate specification" annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="Type of end plate specification",
          __Dymola_joinNext=true));
      Real u_ca_end_plate=333.15*U.K
        "Value of the cathode end plate specification" annotation (Dialog(
          tab="Conditions",
          group="Cathode",
          __Dymola_descriptionLabel=true,
          __Dymola_label="<html>Value (<i>u</i><sub>ca end plate</sub>)</html>"));
      Q.TemperatureAbsolute T_ca "Temperature of cathode end plate";
      Q.Conductance G_ca
        "Thermal conductance of the cathode end plate to the environment";
      Q.Power Qdot_ca "Heat flow rate from the cathode end plate";

      // Material properties
      replaceable package DataH2 = Characteristics.IdealGas constrainedby
        Characteristics.BaseClasses.CharacteristicEOS
        "<html>H<sub>2</sub> gas</html>" annotation (Dialog(tab="Advanced",
            group="Fluid equations of state"), choicesAllMatching=true);
      replaceable package DataH2O = Characteristics.IdealGas constrainedby
        Characteristics.BaseClasses.CharacteristicEOS
        "<html>H<sub>2</sub>O gas</html>" annotation (Dialog(tab="Advanced",
            group="Fluid equations of state"), choicesAllMatching=true);
      replaceable package DataH2Ol = Characteristics.H2O.Liquid constrainedby
        Characteristics.BaseClasses.CharacteristicEOS
        "<html>H<sub>2</sub>O liquid</html>" annotation (Dialog(tab="Advanced",
            group="Fluid equations of state"), choicesAllMatching=true);
      replaceable package DataN2 = Characteristics.IdealGas constrainedby
        Characteristics.BaseClasses.CharacteristicEOS
        "<html>N<sub>2</sub> gas</html>" annotation (Dialog(tab="Advanced",
            group="Fluid equations of state"), choicesAllMatching=true);
      replaceable package DataO2 = Characteristics.IdealGas constrainedby
        Characteristics.BaseClasses.CharacteristicEOS
        "<html>O<sub>2</sub> gas</html>" annotation (Dialog(tab="Advanced",
            group="Fluid equations of state"), choicesAllMatching=true);

      // Standard conditions
      parameter Q.TemperatureAbsolute T_0=273.15*U.K "Temperature" annotation (
          Dialog(
          tab="Advanced",
          group="Standard conditions (for volumetric flow rate)",
          __Dymola_label="<html><i>T</i><sub>0</sub>"));
      parameter Q.PressureAbsolute p_0=U.atm "Pressure" annotation (Dialog(
          tab="Advanced",
          group="Standard conditions (for volumetric flow rate)",
          __Dymola_label="<html><i>p</i><sub>0</sub>"));

      // Derived and measured conditions
      Q.CurrentAreic zJ_seg[n_y, n_z] "Current density of the segments";
      Q.PressureAbsolute p_sat_an_in "Saturation pressure at the anode inlet";
      Q.PressureAbsolute p_sat_ca_in "Saturation pressure at the cathode inlet";
      Q.Current Ndot_H2Ol_an_in "Flow rate of liquid water into anode";
      Q.Current Ndot_H2Ol_ca_in "Flow rate of liquid water into cathode";
      Q.Pressure p_H2Ol_an_in
        "Non-equilibrium pressure on the H2O liquid at the anode inlet";
      Q.Pressure p_H2Ol_ca_in
        "Non-equilibrium pressure on the H2O liquid at the cathode inlet";
      Q.TemperatureAbsolute T_an_out "Anode outlet temperature";
      Q.TemperatureAbsolute T_ca_out "Cathode outlet temperature";
      Q.Velocity phi_an_in[n_y, n_z] "Velocity profile over the anode inlet";
      Q.Velocity phi_ca_in[n_y, n_z] "Velocity profile over the cathode inlet";
      Q.Velocity phi_an_out[n_y, n_z] "Velocity profile over the anode outlet";
      Q.Velocity phi_ca_out[n_y, n_z]
        "Velocity profile over the cathode outlet";

      // Auxiliary measurements
      output Q.Power Wdot(stateSelect=StateSelect.never) = w*zI
        "Electrical power output of the cell";
      output Q.Power Wdot_yz[n_y, n_z](each stateSelect=StateSelect.never) = -
        anBC.graphite.'e-'.face.phi[1] .* anBC.graphite.'e-'.face.mPhidot[1] -
        caBC.graphite.'e-'.face.phi[1] .* caBC.graphite.'e-'.face.mPhidot[1]
        if environment.analysis "Electrical power of the segments";
      output Q.CurrentAreic zJ_yz[n_y, n_z](each stateSelect=StateSelect.never)
         = -anBC.graphite.'e-'.face.phi[1] .* anBC.graphite.'e-'.face.rho if
        environment.analysis "Current densities of the segments";
      output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anSink.gas.H2.face.Ndot
         + anSource.gas.H2.face.Ndot) if environment.analysis
        "Net rate of hydrogen into the cell";
      output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anSink.gas.H2O.face.Ndot
         + anSource.gas.H2O.face.Ndot) + sum(caSink.gas.H2O.face.Ndot +
        caSource.gas.H2O.face.Ndot) if environment.analysis
        "Net rate of water from the cell";
      output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caSink.gas.O2.face.Ndot
         + caSource.gas.O2.face.Ndot) if environment.analysis
        "Net rate of oxygen into the cell";
      output Q.NumberAbsolute anOutletRH(
        stateSelect=StateSelect.never,
        displayUnit="%") = anSink.gas.H2O.materialOut.y/(
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(anSink.gas.H2O.face.T
        /U.K)*U.Pa) if environment.analysis
        "Relative humidity at the anode outlet";
      output Q.NumberAbsolute caOutletRH(
        stateSelect=StateSelect.never,
        displayUnit="%") = caSink.gas.H2O.materialOut.y/(
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(caSink.gas.H2O.face.T
        /U.K)*U.Pa) if environment.analysis
        "Relative humidity at the cathode outlet";

      Connectors.FaceBus an[n_y, n_z] "Interface to the anode end plate"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-120,0}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-160,0})));
      Connectors.FaceBus anNegative[n_x_an, n_z]
        "Negative interface to the anode flow channel" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-80,-100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-40,-160})));
      Connectors.FaceBus anPositive[n_x_an, n_z]
        "Positive interface to the anode flow channel" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-80,100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-40,160})));
      Connectors.FaceBus ca[n_y, n_z] "Interface to the cathode end plate"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={120,0}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={162,0})));
      Connectors.FaceBus caNegative[n_x_ca, n_z]
        "Negative interface to the cathode flow channel" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={80,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={40,-160})));
      Connectors.FaceBus caPositive[n_x_ca, n_z]
        "Positive interface to the cathode flow channel" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={80,100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={40,160})));

      ByConnector.FaceBus.Single.FlowsGraphiteOnly anBC[n_y, n_z](each graphite(
          'incle-'=true,
          'e-'(
            redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare Modelica.Blocks.Sources.Ramp normalSet(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1),
            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=environment.T)),
          'inclC+'=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-104,0})));

      ByConnector.FaceBus.Single.FlowsGraphiteOnly caBC[n_y, n_z](each graphite(
          'incle-'=true,
          'e-'(redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
              thermalSet(y=environment.T)),
          'inclC+'=true)) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={104,0})));
      ByConnector.FaceBus.Single.Efforts anSource[n_x_an, n_z](gas(
          each inclH2=true,
          each inclH2O=true,
          H2(
            redeclare each function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.current,
            materialSet(y=(anSource.gas.H2.normalOut.y - fill(
                      zI*anStoich/(2*A_an),
                      n_x_an,
                      n_z)) .* A_an_seg),
            redeclare each function normalMeas =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare each function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.force,
            each thermalSet(y=environment.T)),
          H2O(
            redeclare each function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.current,
            materialSet(y=(anSource.gas.H2O.normalOut.y - fill(
                      zI*anStoich*environment.p_H2O/(2*(environment.p -
                    environment.p_H2O)*A_an),
                      n_x_an,
                      n_z)) .* A_an_seg),
            redeclare each function normalMeas =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare each function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.force,
            each thermalSet(y=environment.T))),each liquid(inclH2O=true))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-104,-60})));

      ByConnector.FaceBus.Single.Flows anSink[n_x_an, n_z](gas(
          each inclH2=true,
          each inclH2O=true,
          H2(redeclare each function materialMeas =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                  redeclare package Data = FCSys.Characteristics.IdealGas),
              normalSet(y=p_an_out_noneq*A_an_seg .* anSink.gas.H2.face.rho ./
                  (anSink.gas.H2.face.rho + anSink.gas.H2O.face.rho))),
          H2O(redeclare each function materialMeas =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                  redeclare package Data = FCSys.Characteristics.IdealGas),
              normalSet(y=p_an_out_noneq*A_an_seg .* anSink.gas.H2O.face.rho
                   ./ (anSink.gas.H2.face.rho + anSink.gas.H2O.face.rho)))),
          each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-56,60})));
      ByConnector.FaceBus.Single.Efforts caSource[n_x_ca, n_z](gas(
          each inclH2O=true,
          each inclN2=true,
          each inclO2=true,
          H2O(
            redeclare each function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.current,
            materialSet(y=(caSource.gas.H2O.normalOut.y - fill(
                      zI*caStoich*environment.p_H2O/(4*environment.p_O2*A_ca),
                      n_x_ca,
                      n_z)) .* A_ca_seg),
            redeclare each function normalMeas =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare each function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.force,
            each thermalSet(y=environment.T)),
          N2(
            redeclare each function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.current,
            materialSet(y=(caSource.gas.N2.normalOut.y - fill(
                      zI*caStoich*(environment.p - environment.p_O2 -
                    environment.p_H2O)/(4*environment.p_O2*A_ca),
                      n_x_ca,
                      n_z)) .* A_ca_seg),
            redeclare each function normalMeas =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare each function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.force,
            each thermalSet(y=environment.T)),
          O2(
            redeclare each function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.current,
            materialSet(y=(caSource.gas.O2.normalOut.y - fill(
                      zI*caStoich/(4*A_ca),
                      n_x_ca,
                      n_z)) .* A_ca_seg),
            redeclare each function normalMeas =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

            redeclare each function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.force,

            each thermalSet(y=environment.T))),each liquid(inclH2O=true))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={56,-60})));

      ByConnector.FaceBus.Single.Flows caSink[n_x_ca, n_z](gas(
          each inclH2O=true,
          each inclN2=true,
          each inclO2=true,
          H2O(redeclare each function materialMeas =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                  redeclare package Data = FCSys.Characteristics.IdealGas),
              normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.H2O.face.rho
                   ./ (caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho +
                  caSink.gas.O2.face.rho))),
          N2(redeclare each function materialMeas =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                  redeclare package Data = FCSys.Characteristics.IdealGas),
              normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.N2.face.rho ./
                  (caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho + caSink.gas.O2.face.rho))),

          O2(redeclare each function materialMeas =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                  redeclare package Data = FCSys.Characteristics.IdealGas),
              normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.O2.face.rho ./
                  (caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho + caSink.gas.O2.face.rho)))),
          each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={104,60})));

      Modelica.Blocks.Sources.RealExpression anPressSet(y(unit="m/(l.T2)") =
          environment.p) "Setpoint for the total anode outlet pressure"
        annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
      Modelica.Blocks.Sources.RealExpression caPressSet(y(unit="m/(l.T2)") =
          environment.p) "Setpoint for the total cathode outlet pressure"
        annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
      Modelica.Blocks.Sources.RealExpression anThermoPress(y(unit="m/(l.T2)")
           = sum((anSink.gas.H2.materialOut.y + anSink.gas.H2O.materialOut.y)
           .* outerProduct(L_x_an, L_z))/A_an)
        "Thermodynamic pressure at the anode outlet"
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      Modelica.Blocks.Sources.RealExpression caThermoPress(y(unit="m/(l.T2)")
           = sum((caSink.gas.H2O.materialOut.y + caSink.gas.N2.materialOut.y +
          caSink.gas.O2.materialOut.y) .* A_ca_seg)/A_ca)
        "Thermodynamic pressure at the cathode outlet"
        annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));
      Modelica.Blocks.Math.Feedback anNoneqSet
        annotation (Placement(transformation(extent={{-30,10},{-10,30}})));
      Modelica.Blocks.Math.Feedback caNoneqSet
        annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
      Modelica.Blocks.Continuous.FirstOrder anValveDynamics(initType=Modelica.Blocks.Types.Init.InitialOutput,
          T=2) "Dynamics of the anode exit valve"
        annotation (Placement(transformation(extent={{0,10},{20,30}})));
      Modelica.Blocks.Continuous.FirstOrder caValveDynamics(initType=Modelica.Blocks.Types.Init.InitialOutput,
          T=2) "Dynamics of the cathode exit valve"
        annotation (Placement(transformation(extent={{0,-30},{20,-10}})));
    protected
      Connectors.RealOutputInternal p_an_out_noneq(unit="m/(l.T2)")
        "Nonequilibrium pressure at the anode outlet" annotation (Placement(
            transformation(extent={{34,10},{54,30}}), iconTransformation(extent
              ={{174,-40},{194,-20}})));
      Connectors.RealOutputInternal p_ca_out_noneq(unit="m/(l.T2)")
        "Nonequilibrium pressure at the cathode outlet" annotation (Placement(
            transformation(extent={{34,-30},{54,-10}}), iconTransformation(
              extent={{174,-80},{194,-60}})));
    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // **remove liq H2O at inlet and outlet (assume no diffusion, adiabatic, no shear force)

      // Electrical
      w = R*zI;
      P = w*zI;
      A*zJ = zI;
      zJ_seg = anBC.graphite.'e-'.face.phi[1] .* anBC.graphite.'e-'.face.rho;
      zI = sum(zJ_seg .* A_seg);

      // Anode humidity
      p_sat_an_in = saturationPressureSI(T_an_in/U.K)*U.Pa;
      p_H2O_an_in = saturationPressureSI(T_sat_an_in/U.K)*U.Pa;
      p_H2O_an_in = min(anInletRH, 1)*p_sat_an_in;
      // Liquid makes up the remainder if RH > 100%:
      Ndot_H2Ol_an_in = max(anInletRH - 1, 0)*sum(phi_an_in .* A_an_seg)/
        DataH2O.v_Tp(T_an_in, p_sat_an_in);
      Ndot_H2Ol_an_in = sum(anSource.liquid.H2O.face.phi[Orient.normal] .*
        A_an_seg)/DataH2Ol.v_Tp(T_an_in);

      // Cathode humidity
      p_sat_ca_in = saturationPressureSI(T_ca_in/U.K)*U.Pa;
      p_H2O_ca_in = saturationPressureSI(T_sat_ca_in/U.K)*U.Pa;
      p_H2O_ca_in = max(caInletRH, 1)*p_sat_ca_in;
      // Liquid makes up the remainder if RH > 100%:
      Ndot_H2Ol_ca_in = max(caInletRH - 1, 0)*sum(phi_ca_in .* A_ca_seg)/
        DataH2O.v_Tp(T_ca_in, p_sat_ca_in);
      Ndot_H2Ol_ca_in = sum(caSource.liquid.H2O.face.phi[Orient.normal] .*
        A_ca_seg)/DataH2Ol.v_Tp(T_ca_in);

      // End plates
      Qdot_an = G_an*(T_an - environment.T) "Anode";
      Qdot_ca = G_ca*(T_ca - environment.T) "Cathode";
      Qdot_an = sum(anBC.graphite.'C+'.face.Qdot + anBC.graphite.'e-'.face.Qdot);
      Qdot_ca = sum(caBC.graphite.'C+'.face.Qdot + caBC.graphite.'e-'.face.Qdot);

      // Anode flow rate
      anStoich*zI = I_an;
      J_an*A = I_an;
      Vdot_g_an_in = inSign(anInletSide)*sum(outerProduct(L_x_an, L_z) .*
        phi_an_in);
      I_an = sum(phi_an_in .* anSource.gas.H2.face.rho .* A_an_seg);

      // Cathode flow rate
      caStoich*zI = I_ca;
      J_ca*A = I_ca;
      Vdot_g_ca_in = inSign(caInletSide)*sum(outerProduct(L_x_ca, L_z) .*
        phi_ca_in);
      I_ca = sum(phi_ca_in .* caSource.gas.H2O.face.rho .* A_ca_seg);

      // Pressures at the inlets and outlets
      for j in 1:n_z loop
        for i in 1:n_x_an loop
          p_an_in = DataH2.p_Tv(anSource[i, j].gas.H2.face.T, 1/anSource[i, j].gas.H2.face.rho)
             + DataH2O.p_Tv(anSource[i, j].gas.H2O.face.T, 1/anSource[i, j].gas.H2O.face.rho)
             - inSign(anInletSide)*(anSource[i, j].gas.H2.face.mPhidot[Orient.normal]
             + anSource[i, j].gas.H2O.face.mPhidot[Orient.normal])/A_an_seg[i,
            j];
          p_an_out = DataH2.p_Tv(anSink[i, j].gas.H2.face.T, 1/anSink[i, j].gas.H2.face.rho)
             + DataH2O.p_Tv(anSink[i, j].gas.H2O.face.T, 1/anSink[i, j].gas.H2O.face.rho)
             + inSign(anInletSide)*(anSink[i, j].gas.H2.face.mPhidot[Orient.normal]
             + anSink[i, j].gas.H2O.face.mPhidot[Orient.normal])/A_ca_seg[i, j];
        end for;
        for i in 1:n_x_ca loop
          p_ca_in = DataH2O.p_Tv(caSource[i, j].gas.H2O.face.T, 1/caSource[i, j].gas.H2O.face.rho)
             + DataN2.p_Tv(caSource[i, j].gas.N2.face.T, 1/caSource[i, j].gas.N2.face.rho)
             + DataO2.p_Tv(caSource[i, j].gas.O2.face.T, 1/caSource[i, j].gas.O2.face.rho)
             - inSign(caInletSide)*(caSource[i, j].gas.H2O.face.mPhidot[Orient.normal]
             + caSource[i, j].gas.N2.face.mPhidot[Orient.normal] + caSource[i,
            j].gas.O2.face.mPhidot[Orient.normal])/A_ca_seg[i, j];
          p_ca_out = DataH2O.p_Tv(caSink[i, j].gas.H2O.face.T, 1/caSink[i, j].gas.H2O.face.rho)
             + DataN2.p_Tv(caSink[i, j].gas.N2.face.T, 1/caSink[i, j].gas.N2.face.rho)
             + DataO2.p_Tv(caSink[i, j].gas.O2.face.T, 1/caSink[i, j].gas.O2.face.rho)
             + inSign(caInletSide)*(caSink[i, j].gas.H2O.face.mPhidot[Orient.normal]
             + caSink[i, j].gas.N2.face.mPhidot[Orient.normal] + caSink[i, j].gas.O2.face.mPhidot[
            Orient.normal])/A_ca_seg[i, j];
        end for;
      end for;

      // Assumptions
      0 = sum(anSink.gas.H2.face.Qdot + anSink.gas.H2O.face.Qdot + anSink.liquid.H2O.face.Qdot)
        "Adiabatic across the anode outlet";
      0 = sum(caSink.gas.H2O.face.Qdot + caSink.gas.N2.face.Qdot + caSink.gas.O2.face.Qdot
         + caSink.liquid.H2O.face.Qdot) "Adiabatic across the cathode outlet";

      if anInletSide == Side.n then
        connect(anSource.face, anNegative) annotation (Line(
            points={{-100,-60},{-90,-60},{-90,-90},{-80,-100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anSink.face, anPositive) annotation (Line(
            points={{-60,60},{-70,60},{-70,90},{-80,100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
      else
        connect(anSource.face, anPositive) annotation (Line(
            points={{-100,-60},{-90,-60},{-90,90},{-80,100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None,
            pattern=LinePattern.Dash));
        connect(anSink.face, anNegative) annotation (Line(
            points={{-60,60},{-70,60},{-70,-90},{-80,-100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None,
            pattern=LinePattern.Dash));
      end if;
      if caInletSide == Side.n then
        connect(caSource.face, caNegative) annotation (Line(
            points={{60,-60},{70,-60},{70,-90},{80,-100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.face, caPositive) annotation (Line(
            points={{100,60},{90,60},{90,90},{80,100}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
      else
        connect(caSource.face, caPositive) annotation (Line(
            points={{60,-60},{70,-60},{70,90},{80,100}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.face, caNegative) annotation (Line(
            points={{100,60},{90,60},{90,-90},{80,-100}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            thickness=0.5,
            smooth=Smooth.None));
      end if;
      connect(anBC.face, an) annotation (Line(
          points={{-108,1.23436e-15},{-108,5.55112e-16},{-120,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caBC.face, ca) annotation (Line(
          points={{108,-1.34539e-15},{112,0},{114,-7.90278e-16},{114,
              5.55112e-16},{120,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(anValveDynamics.y, p_an_out_noneq) annotation (Line(
          points={{21,20},{44,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(caValveDynamics.y, p_ca_out_noneq) annotation (Line(
          points={{21,-20},{44,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(caPressSet.y, caNoneqSet.u1) annotation (Line(
          points={{-39,-20},{-28,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(caThermoPress.y, caNoneqSet.u2) annotation (Line(
          points={{-39,-40},{-20,-40},{-20,-28}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(caNoneqSet.y, caValveDynamics.u) annotation (Line(
          points={{-11,-20},{-2,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(anPressSet.y, anNoneqSet.u1) annotation (Line(
          points={{-39,20},{-28,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(anNoneqSet.u2, anThermoPress.y) annotation (Line(
          points={{-20,12},{-20,6.10623e-16},{-39,6.10623e-16}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(anNoneqSet.y, anValveDynamics.u) annotation (Line(
          points={{-11,20},{-2,20}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        structurallyIncomplete=true,
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,-100},
                {120,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{
                160,160}}), graphics={Rectangle(
                  extent={{-160,160},{160,-160}},
                  lineColor={191,191,191},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Backward),Rectangle(extent={{-160,160},
              {160,-160}}, lineColor={0,0,0})}),
        Documentation(info="
    <html>
    <p>Any of the settings for the operating conditions can be time-varying expressions.
    In each group,
    specify exactly one variable (otherwise the model will be structurally singular).</p>

    <p>The relative humidity (<code>anInletRH</code> or <code>caInletRH</code>) may specified to be greater 
    than 100 %.  In that case, liquid
    water is injected to provide the amount above saturation.  The relative humidity
    is taken to be equal to the quotient of the H<sub>2</sub>O vapor pressure 
    (<i>p</i><sub>H2O an in</sub> or <i>p</i><sub>H2O ca in</sub>) and the saturation pressure.
    Therefore, liquid water will also be injected if the specified vapor pressure is specified 
    to be above saturation pressure or the specified dew point (<i>T</i><sub>sat an in</sub> or 
    <i>T</i><sub>sat ca in</sub>) is above the actual temperature.</p>
    
    <p><i>Equivalent current</i> is the rate of supply of a reactant required to support the
    given current
    assuming the reactant is entirely consumed (complete utilization).</p>

    <p>Assumptions **review and update:
    <ol>
    <li>The outer x-axis surface of each end plate is each uniform in temperature.</li>
    <li>No heat is conducted from the rest of the cell hardware.</li>
    <li>The voltage is uniform across each end plate.</li>
    <li>There is no turbulence in the fluid at either inlet (i.e., zero transverse velocity
    at each inlet boundary).</li>
    <li>There is no shear force on the fluid at either outlet.</li>
    <li>The species (gases and liquid) of each stream have the same temperature at each inlet and outlet.</li>
    <li>The sum of the thermodynamic and nonequilibrium pressure is uniform over each inlet and outlet.</li>
    <li>The temperature is uniform over each inlet and outlet.</li>
    <li>There is no diffusion of the reactants (H<sub>2</sub> and O<sub>2</sub>) or liquid water
    into the cell (only advection).</li>
    <li>There is no diffusion of the fluid species 
    (H<sub>2</sub>, H<sub>2</sub>O, N<sub>2</sub>, and O<sub>2</sub>) 
    out of the cell (advection only).</li>
    <li>The inlet and outlet pressures are applied to the gas mixture by Dalton's law.</li>
    <li>At the inlet, the liquid has the pressure necessary and sufficient for the prescribed 
    humidity (zero unless RH > 100%).</li>
    <li>At the outlet, the liquid has the same pressure as the gas (Amagat's law).</li>
    <li>There is no net thermal conduction across either outlet.</li>
    </ol></p>
        
    <li>There is no nonequilibrium force on any species at either inlet.  This means that
    the velocity of the first subregion along the channel will be the same is the velocity
    at the inlet.</li>

    <p>The temperatures of the endplates (<i>T</i><sub>an</sub> and <i>T</i><sub>ca</sub>)
    should not be equal to the temperature of the environment unless <i>G</i><sub>an</sub>
    and <i>G</i><sub>ca</sub> are explicitly set.  Otherwise there will be a mathematical
    singularity.  Regard the environment as the ambient conditions, not the conditions to
    which the cell is held.</p>
    </html>"));
    end TestStand;

    package Enumerations "Choices of options"

      extends Modelica.Icons.BasesPackage;

      type ElectricalSpec = enumeration(
          currentDensity "Current density",
          current "Current",
          voltage "Voltage",
          resistance "Resistance",
          power "Power") "Ways to specify the electrical load";
      type FlowSpec = enumeration(
          stoich "Stoichiometric rate",
          currentDensity "Equivalent current density",
          current "Equivalent current",
          volumetric "Standard volumetric rate (conditions on Advanced tab)",
          pressure "Inlet pressure") "Ways to specify the anode flow rate";

      type HumiditySpec = enumeration(
          relative "Relative humidity",
          pressure "Vapor pressure",
          dewPoint "Dew point") "Ways to specify humidity";
      type ThermalSpec = enumeration(
          temperature "Temperature",
          conductance "Thermal conductance with the environment",
          rate "Heat flow rate") "Ways to specify a thermal condition";
    end Enumerations;
  end TestStands;

  record Environment "Environmental properties for a simulation"
    extends FCSys.Icons.Names.Top3;

    // Store the values of the base constants and units.
    final constant U.Bases.Base baseUnits=U.base "Base constants and units"
      annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

    parameter Boolean analysis=true "Include optional variables for analysis"
      annotation (choices(__Dymola_checkBox=true));

    parameter Q.TemperatureAbsolute T(nominal=300*U.K) = 298.15*U.K
      "Temperature" annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));

    parameter Q.PressureAbsolute p(nominal=U.atm) = U.atm "Pressure"
      annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));

    parameter Q.NumberAbsolute RH(
      displayUnit="%",
      max=1) = 0.6 "Relative humidity";
    final parameter Q.PressureAbsolute p_H2O=RH*
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa
      "Pressure of H2O vapor";
    parameter Q.NumberAbsolute n_O2(
      final max=1,
      displayUnit="%") = 0.208
      "<html>Dry-gas concentration of O<sub>2</sub></html>"
      annotation (Dialog(__Dymola_label="<html><i>n</i><sub>O2</sub></html>"));
    // Value from http://en.wikipedia.org/wiki/Oxygen
    final parameter Q.PressureAbsolute p_O2=n_O2*(p - p_H2O) "Pressure of O2";

    parameter Q.Acceleration a[Axis]={0,Modelica.Constants.g_n*U.m/U.s^2,0}
      "Acceleration due to body forces"
      annotation (Dialog(__Dymola_label="<html><b><i>a</i></b></html>"));

    // The gravity component is positive because it's added to the transient
    // term in the Species model.
    parameter Q.ForceSpecific E[Axis]={0,0,0} "Electric field"
      annotation (Dialog(__Dymola_label="<html><b><i>E</i></b></html>"));

    annotation (
      defaultComponentPrefixes="inner",
      missingInnerMessage="Your model is using an outer \"environment\" record, but an inner \"environment\"
record is not defined.  For simulation, drag FCSys.Conditions.Environment into
your model to specify global conditions and defaults.  Otherwise, the default
settings will be used.",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            extent={{-80,60},{80,-100}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-70,50},{70,-98}},
            lineColor={255,255,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={170,213,255}),
          Rectangle(
            extent={{-72,-60},{72,-100}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Line(
            points={{-40,-20},{-30,-28},{-10,-50},{-10,-50},{16,-12},{40,0}},
            color={0,0,0},
            smooth=Smooth.Bezier),
          Ellipse(
            extent={{32,8},{48,-8}},
            pattern=LinePattern.None,
            lineColor={170,213,255},
            fillColor={50,50,50},
            fillPattern=FillPattern.Sphere),
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
            pattern=LinePattern.Dash),
          Line(points={{-70,-60},{70,-60}}, color={0,0,0}),
          Line(points={{-66,-90},{-36,-60}}, color={0,0,0}),
          Line(points={{2,-90},{32,-60}}, color={0,0,0}),
          Line(points={{36,-90},{66,-60}}, color={0,0,0}),
          Line(points={{-32,-90},{-2,-60}}, color={0,0,0})}));

  end Environment;

  model Router "Connect two pairs of faces to pass through or cross over"
    extends FCSys.Icons.Names.Top3;
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
is connected to <code>positive2</code>, as shown by <a href=\"#Fig1a\">Figure 1a</a>.</p>

<p>If <code>crossOver</code> is set to <code>true</code>, then the router will be in cross-over mode.  In that case, <code>negative1</code> is connected to <code>positive2</code>
and <code>negative2</code> is
connected to <code>positive1</code>, as shown by <a href=\"#Fig1b\">Figure 1b</a>.</p>

    <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\" align=center>
      <tr align=center>
        <td id=\"Fig1a\" align=center width=120>
          <img src=\"modelica://FCSys/Resources/Documentation/Conditions/Router/PassThrough.png\">
<br><b>a:</b> Pass-through
        </td>
        <td id=\"Fig1b\" align=center>
          <img src=\"modelica://FCSys/Resources/Documentation/Conditions/Router/CrossOver.png\">
<br><b>b:</b> Cross-over
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
end Conditions;
