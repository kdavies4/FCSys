within FCSys;
package Conditions "Models to specify and measure operating conditions"
  extends Modelica.Icons.SourcesPackage;
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model BoundaryCondition
      "Test the conditions for the boundary of a subregion"
      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      ByConnector.BoundaryBus.Single.Sink boundary
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
      connect(subregion.yPositive, boundary.boundary) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(NumberOfIntervals=5000), Commands(file=
              "Resources/Scripts/Dymola/Conditions.Examples.BoundaryCondition.mos"
            "Conditions.Examples.BoundaryCondition.mos"));
    end BoundaryCondition;

    model BoundaryConditionPhases
      "Test the conditions for the boundary of a subregion with phases"
      import FCSys.Utilities.Coordinates.cartWrap;
      import Modelica.Math.BooleanVectors.countTrue;
      import Modelica.Math.BooleanVectors.enumerate;
      import Modelica.Math.BooleanVectors.index;
      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) =
        ones(3)*U.cm "Length" annotation (Dialog(group="Geometry",
            __Dymola_label="<html><b><i>L</i></b></html>"));
      final inner parameter Q.Volume V=product(L) "Volume";

      // Included boundaries
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

      ByConnector.BoundaryBus.Single.Phases.Gas boundary(inclH2O=true, H2O(
            redeclare
            Conditions.ByConnector.Boundary.Single.ThermalDiffusive.heatRate
            thermal, redeclare
            Conditions.ByConnector.Boundary.Single.Material.Current material(
              set(y=U.A))))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));

      Conditions.ByConnector.Amagat.VolumeFixed volume(n_phases=1)
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      FCSys.Phases.Gas gas(
        inclH2=false,
        inclH2O=true,
        final n_trans=n_trans)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    protected
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ} "true, if each pairs of boundaries is included";
      final inner parameter Boolean inclRot[Axis]={inclTransY and inclTransZ,
          inclTransZ and inclTransX,inclTransX and inclTransY}
        "true, if each axis of rotation has all its tangential boundaries included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartTrans[n_trans]=index(inclTrans)
        "Cartesian-axis indices of the transport axes";
      final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
        "Transport-axis indices of the Cartesian axes";

    equation
      connect(gas.yPositive, boundary.boundary) annotation (Line(
          points={{0,10},{0,16},{6.10623e-016,16},{6.10623e-016,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(volume.dalton[1], gas.dalton) annotation (Line(
          points={{11,-11},{8,-8}},
          color={47,107,251},
          smooth=Smooth.None));
      annotation (experiment);
    end BoundaryConditionPhases;

    model Router
      "<html>Test the <a href=\"modelica://FCSys.Conditions.Router\">Router</a> model</html>"
      extends Modelica.Icons.Example;

      Conditions.Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      ByConnector.BoundaryBus.Single.Source fastFlow(gas(inclH2=true,H2(
              materialSet(y=-U.A)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-30,20})));
      ByConnector.BoundaryBus.Single.Source slowFlow(gas(inclH2=true,H2(
              materialSet(y=-U.mA)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-30,-20})));

      ByConnector.BoundaryBus.Single.Sink sink2(gas(inclH2=true)) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-20})));
      ByConnector.BoundaryBus.Single.Sink sink1(gas(inclH2=true)) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,20})));
    equation
      connect(router.positive2, sink1.boundary) annotation (Line(
          points={{8,4},{20,4},{20,20},{26,20},{26,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(router.positive1, sink2.boundary) annotation (Line(
          points={{8,-4},{20,-4},{20,-20},{26,-20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(fastFlow.boundary, router.negative2) annotation (Line(
          points={{-26,20},{-20,20},{-20,4},{-8,4}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(slowFlow.boundary, router.negative1) annotation (Line(
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
      extends Modelica.Icons.UnderConstruction;
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
          points={{10,-20},{10,2},{-6,2}},
          color={0,0,255},
          smooth=Smooth.None));

      connect(subregion.xPositive, anodeAdapter.boundary) annotation (Line(
          points={{-30,6.10623e-016},{-24,6.10623e-016},{-24,6.10623e-016},{-14,
              6.10623e-016}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(gasVolume.heatPort, anodeAdapter.heatPort) annotation (Line(
          points={{30,30},{20,30},{20,-2},{-6,-2}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(gasVolume.ports[1], anodeAdapter.gasPort) annotation (Line(
          points={{40,20},{40,6},{-6,6}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquidVolume.heatPort, anodeAdapter.heatPort) annotation (Line(
          points={{70,30},{60,30},{60,-2},{-6,-2}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(anodeAdapter.liquidPort, liquidVolume.ports[1]) annotation (Line(
          points={{-6,-6},{80,-6},{80,20}},
          color={0,127,255},
          smooth=Smooth.None));
      annotation (experiment(StopTime=2e-10), Commands(file=
              "Resources/Scripts/Dymola/Conditions.Examples.Adapteminus.mos"
            "Conditions.Examples.Adapteminus.mos"));
    end AnodeAdapter;

    model Stoichiometry "Test the stoichiometry of a reaction"

      extends Modelica.Icons.Example;
      import Modelica.Math.BooleanVectors.countTrue;
      extends Modelica.Icons.UnderConstruction;
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

      replaceable ByConnector.Chemical.Current speciesA(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=2000*U.K,
        redeclare Modelica.Blocks.Sources.Ramp set(duration=100, height=-1*U.A))
        annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
      replaceable ByConnector.Chemical.Potential speciesB(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=3000*U.K)
        annotation (Placement(transformation(extent={{30,0},{50,20}})));
      replaceable ByConnector.Chemical.Potential speciesC(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=4000*U.K)
        annotation (Placement(transformation(extent={{30,-20},{50,0}})));
      Conditions.Adapters.ChemicalReaction A(
        m=U.g/U.mol,
        final n_trans=n_trans,
        n=-1)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      Conditions.Adapters.ChemicalReaction B(
        m=U.g/U.mol,
        final n_trans=n_trans,
        n=2) annotation (Placement(transformation(extent={{30,-10},{10,10}})));

      Conditions.Adapters.ChemicalReaction C(
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
      connect(A.chemical, speciesA.chemical) annotation (Line(
          points={{-24,0},{-40,0},{-40,6}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(B.chemical, speciesB.chemical) annotation (Line(
          points={{24,0},{40,0},{40,6}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(C.chemical, speciesC.chemical) annotation (Line(
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

    model AmagatDalton
      "<html>Adapter between the <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> and <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connectors</html>"

      extends FCSys.Icons.Names.Top1;

      Connectors.Amagat amagat "Connector for additivity of volume" annotation
        (Placement(transformation(extent={{-30,-10},{-10,10}}),
            iconTransformation(extent={{-50,-10},{-30,10}})));
      Connectors.Dalton dalton "Connector for additivity of pressure"
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
              smooth=Smooth.None), Polygon(
              points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
              lineColor={47,107,251},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
    end AmagatDalton;

    model ChemicalReaction
      "<html>Adapter between the <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> and <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connectors</html>"

      // extends FCSys.Icons.Names.Top1;

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
      output Q.Velocity phi[n_trans](each stateSelect=StateSelect.never) =
        actualStream(chemical.phi) if environment.analysis
        "Velocity of the stream";
      output Q.PotentialAbsolute sT(stateSelect=StateSelect.never) =
        actualStream(chemical.sT) if environment.analysis
        "Specific entropy-temperature product of the stream";

      Connectors.Chemical chemical(redeclare final constant Integer n_trans=
            n_trans) "Connector for a species in a chemical reaction"
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}}),
            iconTransformation(extent={{-50,-10},{-30,10}})));
      // Note:  This redeclaration is necessary due to errors in Dymola 2014.
      Connectors.Reaction reaction(final n_trans=n_trans)
        "Connector for a chemical reaction" annotation (Placement(
            transformation(extent={{10,-10},{30,10}}), iconTransformation(
              extent={{30,-10},{50,10}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Equal intensive properties
      reaction.g = n*chemical.g "Chemical potential";
      reaction.phi = chemical.phi "Velocity (upon outflow)";
      reaction.sT = chemical.sT
        "Specific entropy-temperature product (upon outflow)";

      // Conservation (without storage)
      0 = chemical.Ndot + n*reaction.Ndot "Material";
      zeros(n_trans) = m*actualStream(chemical.phi)*chemical.Ndot + reaction.mPhidot
        "Translational momentum";
      0 = actualStream(chemical.sT)*chemical.Ndot + reaction.Qdot "Energy";
      annotation (
        Documentation(info="<html><p>This model is used to add the stoichiometrically-weighted chemical potential
    of a species to the net chemical potential of a reaction.  The species is produced at the
    stoichiometrically-weighted rate of the reaction.</p>

    <p>
    For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={
            Line(
              points={{-30,0},{30,0}},
              color={255,195,38},
              smooth=Smooth.None),
            Text(
              extent={{-100,20},{100,60}},
              lineColor={0,0,0},
              textString="%n %name"),
            Polygon(
              points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
              lineColor={255,195,38},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end ChemicalReaction;

    package MSL
      "<html>Adapters to the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
      extends Modelica.Icons.Package;

      model Anode
        "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the boundary connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>"
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

        Connectors.BoundaryBus boundary
          "Multi-species connector for translational momentum and heat"
          annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
              iconTransformation(extent={{-50,-10},{-30,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
            Medium = GasMedium) "Modelica fluid port for the gas" annotation (
            Placement(transformation(extent={{50,50},{70,70}}),
              iconTransformation(extent={{30,50},{50,70}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{50,10},{70,30}}), iconTransformation(extent={{30,10},{50,30}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{50,
                  -30},{70,-10}}), iconTransformation(extent={{30,-30},{50,-10}})));
        Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final
            package Medium = LiquidMedium) "Modelica fluid port for the liquid"
          annotation (Placement(transformation(extent={{50,-70},{70,-50}}),
              iconTransformation(extent={{30,-70},{50,-50}})));
        Phases.AnodeGas gas(redeclare final package Medium = GasMedium)
          "Gas subadapter"
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));
        Phases.Graphite graphite('inclC+'=true, 'incle-'=true)
          "Graphite subadapter"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
          "Liquid subadapter"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        connect(gas.boundary, boundary.gas) annotation (Line(
            points={{-4,40},{-20,40},{-20,5.55112e-016},{-40,5.55112e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(gasPort, gas.fluidPort) annotation (Line(
            points={{60,60},{30,60},{30,44},{4,44}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(gas.heatPort, heatPort) annotation (Line(
            points={{4,36},{20,36},{20,-20},{60,-20}},
            color={191,0,0},
            smooth=Smooth.None));

        connect(graphite.boundary, boundary.graphite) annotation (Line(
            points={{-4,0},{-20,0},{-20,5.55112e-016},{-40,5.55112e-016}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5));

        connect(graphite.pin, pin) annotation (Line(
            points={{4,2},{30,2},{30,20},{60,20}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(graphite.heatPort, heatPort) annotation (Line(
            points={{4,-2},{20,-2},{20,-20},{60,-20}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(liquid.boundary, boundary.liquid) annotation (Line(
            points={{-4,-40},{-20,-40},{-20,5.55112e-016},{-40,5.55112e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(liquidPort, liquid.fluidPort) annotation (Line(
            points={{60,-60},{30,-60},{30,-36},{4,-36}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(liquid.heatPort, heatPort) annotation (Line(
            points={{4,-44},{20,-44},{20,-20},{60,-20}},
            color={191,0,0},
            smooth=Smooth.None));

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={
              Line(
                points={{-40,0},{-20,0}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{0,60},{40,60}},
                color={0,127,255},
                smooth=Smooth.None),
              Line(
                points={{-20,60},{-20,-60}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{0,-20},{40,-20}},
                color={140,0,0},
                smooth=Smooth.None),
              Polygon(
                points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                lineColor={140,0,0},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{0,80},{-20,60},{0,40},{20,60},{0,80}},
                lineColor={0,127,255},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{0,-60},{40,-60}},
                color={0,127,255},
                smooth=Smooth.None),
              Polygon(
                points={{0,-40},{-20,-60},{0,-80},{20,-60},{0,-40}},
                lineColor={0,127,255},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{10,20},{40,20}},
                color={0,0,255},
                visible='incle-',
                smooth=Smooth.None),
              Polygon(
                points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                lineColor={0,0,255},
                visible='incle-',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-20,-48},{20,-72}},
                lineColor={127,127,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                textString="l"),
              Text(
                extent={{-20,72},{20,48}},
                lineColor={127,127,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                textString="g")}), Diagram(coordinateSystem(preserveAspectRatio
                =false, extent={{-100,-100},{100,100}}), graphics));
      end Anode;

      model Cathode
        "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the boundary connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>"
        extends FCSys.Icons.Names.Top4;

        replaceable package GasMedium = Media.CathodeGas constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the gas"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));
        replaceable package LiquidMedium =
            Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Connectors.BoundaryBus boundary
          "Multi-species connector for translational momentum and heat"
          annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
              iconTransformation(extent={{-50,-10},{-30,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
            Medium = GasMedium) "Modelica fluid port for the gas" annotation (
            Placement(transformation(extent={{50,70},{70,90}}),
              iconTransformation(extent={{30,50},{50,70}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{50,30},{70,50}}), iconTransformation(extent={{30,10},{50,30}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
          "Modelica heat port" annotation (Placement(transformation(extent={{50,
                  -12},{70,8}}), iconTransformation(extent={{30,-30},{50,-10}})));
        Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final
            package Medium = LiquidMedium) "Modelica fluid port for the liquid"
          annotation (Placement(transformation(extent={{50,-46},{70,-26}}),
              iconTransformation(extent={{30,-70},{50,-50}})));

        Phases.CathodeGas gas(redeclare final package Medium = GasMedium)
          "Gas subadapter"
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));
        Phases.Graphite graphite('inclC+'=true, 'incle-'=true)
          "Graphite subadapter"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
          "Liquid subadapter"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        connect(gas.boundary, boundary.gas) annotation (Line(
            points={{-4,40},{-20,40},{-20,5.55112e-016},{-40,5.55112e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(gasPort, gas.fluidPort) annotation (Line(
            points={{60,80},{30,80},{30,44},{4,44}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(gas.heatPort, heatPort) annotation (Line(
            points={{4,36},{20,36},{20,-2},{60,-2}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(graphite.boundary, boundary.graphite) annotation (Line(
            points={{-4,6.10623e-016},{-20,6.10623e-016},{-20,5.55112e-016},{-40,
                5.55112e-016}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5));

        connect(graphite.pin, pin) annotation (Line(
            points={{4,2},{40,2},{40,40},{60,40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(graphite.heatPort, heatPort) annotation (Line(
            points={{4,-2},{60,-2}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(liquid.boundary, boundary.liquid) annotation (Line(
            points={{-4,-40},{-20,-40},{-20,5.55112e-016},{-40,5.55112e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(liquidPort, liquid.fluidPort) annotation (Line(
            points={{60,-36},{4,-36}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(liquid.heatPort, heatPort) annotation (Line(
            points={{4,-44},{20,-44},{20,-2},{60,-2}},
            color={191,0,0},
            smooth=Smooth.None));

        annotation (Icon(graphics={
              Line(
                points={{-40,0},{-20,0}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{0,60},{40,60}},
                color={0,127,255},
                smooth=Smooth.None),
              Line(
                points={{-20,60},{-20,-60}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{0,-20},{40,-20}},
                color={140,0,0},
                smooth=Smooth.None),
              Polygon(
                points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                lineColor={140,0,0},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{0,80},{-20,60},{0,40},{20,60},{0,80}},
                lineColor={0,127,255},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{0,-60},{40,-60}},
                color={0,127,255},
                smooth=Smooth.None),
              Polygon(
                points={{0,-40},{-20,-60},{0,-80},{20,-60},{0,-40}},
                lineColor={0,127,255},
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{10,20},{40,20}},
                color={0,0,255},
                visible='incle-',
                smooth=Smooth.None),
              Polygon(
                points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                lineColor={0,0,255},
                visible='incle-',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-20,-48},{20,-72}},
                lineColor={127,127,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                textString="l"),
              Text(
                extent={{-20,72},{20,48}},
                lineColor={127,127,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                textString="g")}), Diagram(coordinateSystem(preserveAspectRatio
                =false, extent={{-100,-100},{100,100}}), graphics));
      end Cathode;

      model Conductor
        "<html>Adapter for electrical and thermal conduction between <a href=\"modelica://Modelica\">Modelica</a> and <a href=\"modelica://FCSys\">FCSys</a></html>"
        extends FCSys.Icons.Names.Top2;

        parameter Boolean 'inclC+'=false "Include C+" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>"));
        parameter Boolean 'incle-'=false "Include e-" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>"));

        Connectors.BoundaryBus boundary "Multi-species connector" annotation (
            Placement(transformation(extent={{-50,-10},{-30,10}}),
              iconTransformation(extent={{-50,-10},{-30,10}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin if 'incle-'
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{30,10},{50,30}}), iconTransformation(extent={{30,10},{50,30}})));
        Phases.Graphite graphite(final 'inclC+'='inclC+', final 'incle-'=
              'incle-') "Graphite subadapter"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort if
          'inclC+' "Modelica heat port" annotation (Placement(transformation(
                extent={{30,-30},{50,-10}}), iconTransformation(extent={{30,-30},
                  {50,-10}})));
      equation
        connect(graphite.boundary, boundary.graphite) annotation (Line(
            points={{-4,6.10623e-016},{-40,6.10623e-016},{-40,0}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5));

        connect(graphite.pin, pin) annotation (Line(
            points={{4,2},{4,2},{20,2},{20,20},{40,20}},
            color={0,0,255},
            smooth=Smooth.None));

        connect(graphite.heatPort, heatPort) annotation (Line(
            points={{4,-2},{20,-2},{20,-20},{40,-20}},
            color={191,0,0},
            smooth=Smooth.None));
        annotation (Icon(graphics={
              Line(
                points={{0,20},{30,20}},
                color={0,0,255},
                visible='incle-',
                smooth=Smooth.None),
              Line(
                points={{-20,20},{-20,0}},
                color={127,127,127},
                visible='incle-',
                smooth=Smooth.None),
              Line(
                points={{-20,0},{-20,-20}},
                color={127,127,127},
                visible='inclC+',
                smooth=Smooth.None),
              Line(
                points={{0,-20},{30,-20}},
                color={140,0,0},
                visible='inclC+',
                smooth=Smooth.None),
              Polygon(
                points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                lineColor={0,0,255},
                visible='incle-',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                lineColor={140,0,0},
                visible='inclC+',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-30,0},{-20,0}},
                color={127,127,127},
                visible='incle-' or 'inclC+',
                smooth=Smooth.None),
              Line(
                points={{0,-20},{40,-20}},
                color={140,0,0},
                visible='inclC+',
                smooth=Smooth.None),
              Line(
                points={{10,20},{40,20}},
                color={0,0,255},
                visible='incle-',
                smooth=Smooth.None),
              Polygon(
                points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                lineColor={140,0,0},
                visible='inclC+',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                lineColor={0,0,255},
                visible='incle-',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-20,0},{-20,20}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{-20,-20},{-20,0}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5),
              Line(
                points={{-40,0},{-20,0}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5)}));
      end Conductor;

      model Electronic
        "<html>Adapter for e<sup>-</sup> between <a href=\"modelica://Modelica\">Modelica</a> and <a href=\"modelica://FCSys\">FCSys</a></html>"
        extends FCSys.Icons.Names.Top1;

        Connectors.BoundaryBus boundary "Multi-species connector" annotation (
            Placement(transformation(extent={{-50,-10},{-30,10}}),
              iconTransformation(extent={{-50,-10},{-30,10}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{30,-10},{50,10}}), iconTransformation(extent={{30,-10},{50,
                  10}})));
        Phases.Graphite graphite(final 'incle-'=true, final 'inclC+'=false)
          "Graphite subadapter"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        connect(graphite.boundary, boundary.graphite) annotation (Line(
            points={{-4,6.10623e-016},{-40,6.10623e-016},{-40,0}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5));

        connect(graphite.pin, pin) annotation (Line(
            points={{4,2},{4,2},{4,0},{40,0}},
            color={0,0,255},
            smooth=Smooth.None));

        annotation (Icon(graphics={
              Line(
                points={{-30,0},{-20,0}},
                color={127,127,127},
                visible='incle-' or 'inclC+',
                smooth=Smooth.None),
              Line(
                points={{10,0},{40,0}},
                color={0,0,255},
                visible='incle-',
                smooth=Smooth.None),
              Polygon(
                points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                lineColor={0,0,255},
                visible='incle-',
                smooth=Smooth.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-40,0},{-20,0}},
                color={127,127,127},
                smooth=Smooth.None,
                thickness=0.5)}), Diagram(coordinateSystem(preserveAspectRatio=
                  false, extent={{-100,-100},{100,100}}), graphics));
      end Electronic;

      package Phases "Adapters for material phases"
        extends Modelica.Icons.Package;

        model AnodeGas
          "<html>Adapter for PEMFC anode gas between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.Icons.Names.Top3;

          replaceable package Medium = Media.AnodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Domains.FluidNeutral H2(redeclare package Medium =
                Modelica.Media.IdealGases.SingleGases.H2, redeclare package
              Data = Characteristics.H2.Gas)
            annotation (Placement(transformation(extent={{-10,10},{10,30}})));
          Domains.FluidNeutral H2O(redeclare package Data =
                Characteristics.H2O.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.H2O)
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
          Junctions.Junction2 junction
            annotation (Placement(transformation(extent={{46,30},{26,50}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{50,30},{70,50}}),
                iconTransformation(extent={{30,30},{50,50}})));
          Connectors.BoundaryBus boundary "FCSys boundary connector"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    50,-50},{70,-30}}), iconTransformation(extent={{30,-50},{50,
                    -30}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Orient]
            "Modelica translational flanges for shear force" annotation (
              Placement(transformation(extent={{50,-10},{70,10}}),
                iconTransformation(extent={{30,-10},{50,10}})));

        equation
          // H2
          connect(H2.boundary, boundary.H2) annotation (Line(
              points={{-4,20},{-20,20},{-20,5.55112e-016},{-40,5.55112e-016}},
              color={127,127,127},
              smooth=Smooth.None));

          connect(H2.heatPort, heatPort) annotation (Line(
              points={{4,16},{40,16},{40,-40},{60,-40}},
              color={191,0,0},
              smooth=Smooth.None));
          connect(H2.fluidPort, junction.purePort1) annotation (Line(
              points={{4,24},{16,24},{16,44},{32,44}},
              color={0,127,255},
              smooth=Smooth.None));

          // H2O
          connect(H2O.boundary, boundary.H2O) annotation (Line(
              points={{-4,-20},{-20,-20},{-20,5.55112e-016},{-40,5.55112e-016}},

              color={127,127,127},
              smooth=Smooth.None));

          connect(H2O.heatPort, heatPort) annotation (Line(
              points={{4,-24},{40,-24},{40,-40},{60,-40}},
              color={191,0,0},
              smooth=Smooth.None));
          connect(H2O.fluidPort, junction.purePort2) annotation (Line(
              points={{4,-16},{20,-16},{20,36},{32,36}},
              color={0,127,255},
              smooth=Smooth.None));

          // Mixture
          connect(junction.mixturePort, fluidPort) annotation (Line(
              points={{40,40},{60,40}},
              color={0,127,255},
              smooth=Smooth.None));
          connect(H2.flange, flange) annotation (Line(
              points={{4,20},{30,20},{30,5.55112e-016},{60,5.55112e-016}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(H2O.flange, flange) annotation (Line(
              points={{4,-20},{30,-20},{30,5.55112e-016},{60,5.55112e-016}},
              color={0,127,0},
              smooth=Smooth.None));
          annotation (Icon(graphics={Line(
                          points={{-40,0},{0,0}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,40},{40,40}},
                          color={0,127,255},
                          smooth=Smooth.None),Line(
                          points={{-20,40},{-20,-40}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,0},{40,0}},
                          color={0,127,0},
                          smooth=Smooth.None),Line(
                          points={{0,-40},{40,-40}},
                          color={140,0,0},
                          smooth=Smooth.None),Polygon(
                          points={{0,-20},{-20,-40},{0,-60},{20,-40},{0,-20}},
                          lineColor={140,0,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                          lineColor={0,127,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,60},{-20,40},{0,20},{20,40},{0,60}},
                          lineColor={0,127,255},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid)}), Diagram(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end AnodeGas;

        model CathodeGas
          "<html>Adapter for PEMFC cathode gas between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.Icons.Names.Top3;

          replaceable package Medium = Media.CathodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Junctions.Junction3 junction(
            redeclare package Medium1 =
                Modelica.Media.IdealGases.SingleGases.H2O,
            redeclare package Medium2 =
                Modelica.Media.IdealGases.SingleGases.N2,
            redeclare package Medium3 =
                Modelica.Media.IdealGases.SingleGases.O2,
            redeclare package MixtureMedium = Medium)
            annotation (Placement(transformation(extent={{46,30},{26,50}})));

          Domains.FluidNeutral H2O(redeclare package Data =
                Characteristics.H2O.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.H2O)
            annotation (Placement(transformation(extent={{-10,10},{10,30}})));
          Domains.FluidNeutral N2(redeclare package Data =
                Characteristics.N2.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.N2)
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
          Domains.FluidNeutral O2(redeclare package Data =
                Characteristics.O2.Gas, redeclare final package Medium =
                Modelica.Media.IdealGases.SingleGases.O2)
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{50,30},{70,50}}),
                iconTransformation(extent={{30,30},{50,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Orient]
            "Modelica translational flanges" annotation (Placement(
                transformation(extent={{50,-10},{70,10}}), iconTransformation(
                  extent={{30,-10},{50,10}})));
          Connectors.BoundaryBus boundary
            "FCSys boundary connector for shear force" annotation (Placement(
                transformation(extent={{-50,-10},{-30,10}}), iconTransformation(
                  extent={{-50,-10},{-30,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    50,-50},{70,-30}}), iconTransformation(extent={{30,-50},{50,
                    -30}})));

        equation
          // H2O
          connect(H2O.boundary, boundary.H2O) annotation (Line(
              points={{-4,20},{-20,20},{-20,5.55112e-016},{-40,5.55112e-016}},
              color={127,127,127},
              smooth=Smooth.None));

          connect(H2O.fluidPort, junction.purePort1) annotation (Line(
              points={{4,24},{12,24},{12,44},{32,44}},
              color={0,127,255},
              smooth=Smooth.None));
          connect(H2O.heatPort, heatPort) annotation (Line(
              points={{4,16},{40,16},{40,-40},{60,-40}},
              color={191,0,0},
              smooth=Smooth.None));

          // N2
          connect(N2.boundary, boundary.N2) annotation (Line(
              points={{-4,6.10623e-016},{-20,6.10623e-016},{-20,5.55112e-016},{
                  -40,5.55112e-016}},
              color={127,127,127},
              smooth=Smooth.None));

          connect(N2.fluidPort, junction.purePort2) annotation (Line(
              points={{4,4},{16,4},{16,40},{32,40}},
              color={0,127,255},
              smooth=Smooth.None));

          connect(N2.heatPort, heatPort) annotation (Line(
              points={{4,-4},{40,-4},{40,-40},{60,-40}},
              color={191,0,0},
              smooth=Smooth.None));

          // O2
          connect(O2.boundary, boundary.O2) annotation (Line(
              points={{-4,-20},{-20,-20},{-20,5.55112e-016},{-40,5.55112e-016}},

              color={127,127,127},
              smooth=Smooth.None));

          connect(O2.fluidPort, junction.purePort3) annotation (Line(
              points={{4,-16},{20,-16},{20,36},{32,36}},
              color={0,127,255},
              smooth=Smooth.None));
          connect(O2.heatPort, heatPort) annotation (Line(
              points={{4,-24},{40,-24},{40,-40},{60,-40}},
              color={191,0,0},
              smooth=Smooth.None));

          // Mixture
          connect(junction.mixturePort, fluidPort) annotation (Line(
              points={{40,40},{60,40}},
              color={0,127,255},
              smooth=Smooth.None));
          connect(H2O.flange, flange) annotation (Line(
              points={{4,20},{30,20},{30,5.55112e-016},{60,5.55112e-016}},
              color={0,127,0},
              smooth=Smooth.None));
          connect(N2.flange, flange) annotation (Line(
              points={{4,6.10623e-016},{24,6.10623e-016},{24,5.55112e-016},{60,
                  5.55112e-016}},
              color={0,127,0},
              smooth=Smooth.None));

          connect(O2.flange, flange) annotation (Line(
              points={{4,-20},{30,-20},{30,5.55112e-016},{60,5.55112e-016}},
              color={0,127,0},
              smooth=Smooth.None));
          annotation (Icon(graphics={Line(
                          points={{-40,0},{0,0}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,40},{40,40}},
                          color={0,127,255},
                          smooth=Smooth.None),Line(
                          points={{-20,40},{-20,-40}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,0},{40,0}},
                          color={0,127,0},
                          smooth=Smooth.None),Line(
                          points={{0,-40},{40,-40}},
                          color={140,0,0},
                          smooth=Smooth.None),Polygon(
                          points={{0,-20},{-20,-40},{0,-60},{20,-40},{0,-20}},
                          lineColor={140,0,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                          lineColor={0,127,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,60},{-20,40},{0,20},{20,40},{0,60}},
                          lineColor={0,127,255},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid)}), Diagram(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end CathodeGas;

        model Graphite
          "<html>Adapter for graphite between <a href=\"modelica://Modelica\">Modelica</a> and <a href=\"modelica://FCSys\">FCSys</a></html>"
          extends FCSys.Icons.Names.Top2;

          parameter Boolean 'inclC+'=false "Include C+" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              group="Species",
              __Dymola_descriptionLabel=false,
              __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>"));
          Domains.Thermal 'C+' if 'inclC+'
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
          parameter Boolean 'incle-'=true "Include e-" annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              group="Species",
              __Dymola_descriptionLabel=true,
              __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>"));
          Domains.Electrical 'e-'(redeclare package Data =
                Characteristics.'e-'.Graphite) if 'incle-'
            annotation (Placement(transformation(extent={{-10,10},{10,30}})));

          Connectors.BoundaryBus boundary "FCSys boundary connector"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort if
            'inclC+' "Modelica heat port" annotation (Placement(transformation(
                  extent={{30,-30},{50,-10}}), iconTransformation(extent={{30,-30},
                    {50,-10}})));
          Modelica.Electrical.Analog.Interfaces.NegativePin pin if 'incle-'
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{30,10},{50,30}}), iconTransformation(extent={{30,10},
                    {50,30}})));

        equation
          // C
          connect('C+'.boundary, boundary.'C+') annotation (Line(
              points={{-4,-20},{-20,-20},{-20,0},{-40,0},{-40,5.55112e-016}},
              color={127,127,127},
              smooth=Smooth.None));

          connect('C+'.heatPort, heatPort) annotation (Line(
              points={{4,-20},{40,-20},{40,-20},{40,-20}},
              color={191,0,0},
              smooth=Smooth.None));

          // e-
          connect('e-'.boundary, boundary.'e-') annotation (Line(
              points={{-4,20},{-20,20},{-20,0},{-40,0},{-40,0}},
              color={127,127,127},
              smooth=Smooth.None));

          connect('e-'.pin, pin) annotation (Line(
              points={{4,20},{40,20}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Icon(graphics={Line(
                          points={{0,20},{30,20}},
                          color={0,0,255},
                          visible='incle-',
                          smooth=Smooth.None),Line(
                          points={{-20,20},{-20,0}},
                          color={127,127,127},
                          visible='incle-',
                          smooth=Smooth.None),Line(
                          points={{-20,0},{-20,-20}},
                          color={127,127,127},
                          visible='inclC+',
                          smooth=Smooth.None),Line(
                          points={{0,-20},{30,-20}},
                          color={140,0,0},
                          visible='inclC+',
                          smooth=Smooth.None),Polygon(
                          points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                          lineColor={0,0,255},
                          visible='incle-',
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                          lineColor={140,0,0},
                          visible='inclC+',
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Line(
                          points={{-30,0},{-20,0}},
                          color={127,127,127},
                          visible='incle-' or 'inclC+',
                          smooth=Smooth.None),Line(
                          points={{0,-20},{40,-20}},
                          color={140,0,0},
                          visible='inclC+',
                          smooth=Smooth.None),Line(
                          points={{10,20},{40,20}},
                          color={0,0,255},
                          visible='incle-',
                          smooth=Smooth.None),Polygon(
                          points={{0,0},{-20,-20},{0,-40},{20,-20},{0,0}},
                          lineColor={140,0,0},
                          visible='inclC+',
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,40},{-20,20},{0,0},{20,20},{0,40}},
                          lineColor={0,0,255},
                          visible='incle-',
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Line(
                          points={{-20,0},{-20,20}},
                          color={127,127,127},
                          visible='incle-',
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{-20,-20},{-20,0}},
                          color={127,127,127},
                          visible='inclC+',
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{-40,0},{-20,0}},
                          color={127,127,127},
                          visible='inclC+' or 'incle-',
                          smooth=Smooth.None,
                          thickness=0.5)}), Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Graphite;

        model Liquid
          "<html>Adapter for liquid between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.Icons.Names.Top3;

          replaceable package Medium =
              Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "Medium model (Modelica)" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Domains.FluidNeutral H2O(redeclare package Data =
                Characteristics.H2O.Liquid, redeclare final package Medium =
                Medium)
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{30,30},{50,50}}),
                iconTransformation(extent={{30,30},{50,50}})));

          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Orient]
            "Modelica translational flanges for shear force" annotation (
              Placement(transformation(extent={{30,-10},{50,10}}),
                iconTransformation(extent={{30,-10},{50,10}})));
          Connectors.BoundaryBus boundary "FCSys boundary connector"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    30,-50},{50,-30}}), iconTransformation(extent={{30,-50},{50,
                    -30}})));

        equation
          // H2O
          connect(H2O.boundary, boundary.H2) annotation (Line(
              points={{-4,6.10623e-016},{-18,6.10623e-016},{-18,0},{-30,0},{-30,
                  5.55112e-016},{-40,5.55112e-016}},
              color={0,0,0},
              smooth=Smooth.None));

          connect(H2O.heatPort, heatPort) annotation (Line(
              points={{4,-4},{20,-4},{20,-40},{40,-40}},
              color={191,0,0},
              smooth=Smooth.None));

          connect(H2O.fluidPort, fluidPort) annotation (Line(
              points={{4,4},{20,4},{20,40},{40,40}},
              color={0,127,255},
              smooth=Smooth.None));
          connect(flange, H2O.flange) annotation (Line(
              points={{40,5.55112e-016},{24,5.55112e-016},{24,6.10623e-016},{4,
                  6.10623e-016}},
              color={0,127,0},
              smooth=Smooth.None));

          annotation (Icon(graphics={Line(
                          points={{-40,0},{0,0}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,40},{40,40}},
                          color={0,127,255},
                          smooth=Smooth.None),Line(
                          points={{-20,40},{-20,-40}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,0},{40,0}},
                          color={0,127,0},
                          smooth=Smooth.None),Line(
                          points={{0,-40},{40,-40}},
                          color={140,0,0},
                          smooth=Smooth.None),Polygon(
                          points={{0,-20},{-20,-40},{0,-60},{20,-40},{0,-20}},
                          lineColor={140,0,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                          lineColor={0,127,0},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid),Polygon(
                          points={{0,60},{-20,40},{0,20},{20,40},{0,60}},
                          lineColor={0,127,255},
                          smooth=Smooth.None,
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid)}), Diagram(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end Liquid;

      end Phases;

      package Domains "Adapters for physical domains"
        extends Modelica.Icons.Package;

        model Electrical
          "<html>Adapter between <a href=\"modelica://Modelica.Electrical.Analog\">Modelica.Electrical.Analog</a> and <a href=\"modelica://FCSys\">FCSys</a></html>"
          import assert = FCSys.Utilities.assertEval;
          extends FCSys.Icons.Names.Top1;

          replaceable package Data = Characteristics.'e-'.Graphite
            constrainedby Characteristics.BaseClasses.Characteristic
            "Characteristic data (for FCSys)" annotation (
            Dialog(group="Material properties"),
            __Dymola_choicesAllMatching=true,
            Placement(transformation(extent={{-60,40},{-40,60}}),
                iconTransformation(extent={{-10,90},{10,110}})));

          Connectors.Boundary boundary "Interface to electrical species"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Electrical.Analog.Interfaces.NegativePin pin
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{30,-10},{50,10}}), iconTransformation(extent={{30,-10},
                    {50,10}})));

        initial equation
          assert(Data.z <> 0, "The species must have charge.");

        equation
          // Assumptions
          boundary.mPhidot = {0,0} "No shear force (assumption #1)";
          boundary.Qdot = 0 "No thermal conduction (assumption #2)";

          Data.g(boundary.T, boundary.p) = Data.z*pin.v*U.V
            "Equal potentials (also conservation of energy)";
          0 = boundary.Ndot + pin.i*U.A/Data.z
            "Conservation of material (also charge), without storage";
          annotation (
            Documentation(info="<html><p>Assumptions:</p><ol>
  <li>There is no shear force  across the interface.</li>
   <li>There is no thermal conduction across the interface.</li>
  </ol>
  <p>Note that the same assumptions are applied in <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a>.
  </p></html>"),
            Icon(graphics={
                Line(
                  points={{10,0},{40,0}},
                  color={0,0,255},
                  smooth=Smooth.None),
                Line(
                  points={{-40,0},{0,0}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={0,0,255},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end Electrical;

        model Thermal
          "<html>Adapter between <a href=\"modelica://Modelica.Thermal.HeatTransfer\">Modelica.Thermal.HeatTransfer</a> and <a href=\"modelica://FCSys\">FCSys</a></html>"
          extends FCSys.Icons.Names.Top1;

          Connectors.ThermalDiffusive boundary
            "Connector for thermal diffusion" annotation (Placement(
                transformation(extent={{-50,-10},{-30,10}}), iconTransformation(
                  extent={{-50,-10},{-30,10}})));

          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    30,-10},{50,10}}), iconTransformation(extent={{30,-10},{50,
                    10}})));

        equation
          boundary.T = heatPort.T*U.K "Equal temperatures";
          0 = boundary.Qdot + heatPort.Q_flow*U.W
            "Conservation of energy, without storage";

          annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={
                    {-100,-100},{100,100}}), graphics={
                Line(
                  points={{-40,0},{-10,0}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Line(
                  points={{10,0},{40,0}},
                  color={140,0,0},
                  smooth=Smooth.None),
                Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={140,0,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Thermal;

        model Fluid
          "<html>Adapter to connect a single fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.Icons.Names.Top3;

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

          Connectors.Boundary boundary
            "Connector for material, momentum, and energy of a single species"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{30,30},{50,50}}),
                iconTransformation(extent={{30,30},{50,50}})));
          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Orient]
            "Modelica translational flanges for shear force" annotation (
              Placement(transformation(extent={{30,-10},{50,10}}),
                iconTransformation(extent={{30,-10},{50,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    30,-50},{50,-30}}), iconTransformation(extent={{30,-50},{50,
                    -30}})));
          Modelica.Electrical.Analog.Interfaces.NegativePin pin
            "Modelica electrical pin" annotation (Placement(transformation(
                  extent={{30,-90},{50,-70}}), iconTransformation(extent={{30,-90},
                    {50,-70}})));

        equation
          // Equal properties
          boundary.p = fluidPort.p*U.Pa "Pressure";
          boundary.phi = der(flange.s)*U.m/U.s "Velocity";
          boundary.T = heatPort.T*U.K "Temperature";
          Medium.specificEnthalpy_pT(fluidPort.p, heatPort.T) = fluidPort.h_outflow;

          // Conservation (without storage)
          0 = Data.z*boundary.Ndot + pin.i*U.A "Charge";
          0 = boundary.Ndot + (fluidPort.m_flow/Data.m)*U.kg/U.s "Material";
          {0,0} = boundary.mPhidot + flange.f*U.N
            "Transverse translational momentum";
          0 = boundary.Qdot + heatPort.Q_flow*U.W "Energy";
          // Note:  All of the advective terms (for all the balance equations)
          // cancel across the interface.
          annotation (
            Documentation(info="<html><p>Assumptions:</p><ol>
  <li>There is no shear force across the interface.</li>
  </ol></html>"),
            Icon(graphics={
                Line(
                  points={{0,40},{40,40}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Line(
                  points={{-20,40},{-20,-80}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Line(
                  points={{0,0},{40,0}},
                  color={0,127,0},
                  smooth=Smooth.None),
                Line(
                  points={{0,-80},{40,-80}},
                  color={0,0,255},
                  smooth=Smooth.None),
                Line(
                  points={{0,-40},{40,-40}},
                  color={140,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{-40,0},{0,0}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Polygon(
                  points={{0,-20},{-20,-40},{0,-60},{20,-40},{0,-20}},
                  lineColor={140,0,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{0,-60},{-20,-80},{0,-100},{20,-80},{0,-60}},
                  lineColor={0,0,255},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={0,127,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{0,60},{-20,40},{0,20},{20,40},{0,60}},
                  lineColor={0,127,255},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end Fluid;

        model FluidNeutral
          "<html>Adapter to connect a single neutral fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          import assert = FCSys.Utilities.assertEval;
          extends FCSys.Icons.Names.Top3;

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

          Connectors.Boundary boundary
            "Connector for material, momentum, and energy of a single species"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final
              package Medium = Medium) "Modelica fluid port" annotation (
              Placement(transformation(extent={{30,30},{50,50}}),
                iconTransformation(extent={{30,30},{50,50}})));
          Modelica.Mechanics.Translational.Interfaces.Flange_a flange[Orient]
            "Modelica translational flanges for shear force" annotation (
              Placement(transformation(extent={{30,-10},{50,10}}),
                iconTransformation(extent={{30,-10},{50,10}})));
          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    30,-50},{50,-30}}), iconTransformation(extent={{30,-50},{50,
                    -30}})));

        initial equation
          assert(Data.z == 0,
            "The species must be neutral, but its chemical formula is " + Data.formula);

        equation
          // Equal properties
          boundary.p = fluidPort.p*U.Pa "Pressure";
          boundary.phi = der(flange.s)*U.m/U.s "Velocity";
          boundary.T = heatPort.T*U.K "Temperature";
          Medium.specificEnthalpy_pT(fluidPort.p, heatPort.T) = fluidPort.h_outflow;

          // Conservation (without storage)
          0 = boundary.Ndot + (fluidPort.m_flow/Data.m)*U.kg/U.s "Material";
          {0,0} = boundary.mPhidot + flange.f*U.N "Translational momentum";
          0 = boundary.Qdot + heatPort.Q_flow*U.W "Energy";
          // Note:  All of the advective terms (for all the balance equations)
          // cancel across the interface.
          annotation (Icon(graphics={
                Line(
                  points={{0,40},{40,40}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Line(
                  points={{0,0},{40,0}},
                  color={0,127,0},
                  smooth=Smooth.None),
                Line(
                  points={{0,-40},{40,-40}},
                  color={140,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{-20,40},{-20,-40}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Line(
                  points={{-40,0},{0,0}},
                  color={127,127,127},
                  smooth=Smooth.None),
                Polygon(
                  points={{0,-20},{-20,-40},{0,-60},{20,-40},{0,-20}},
                  lineColor={140,0,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{0,20},{-20,0},{0,-20},{20,0},{0,20}},
                  lineColor={0,127,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{0,60},{-20,40},{0,20},{20,40},{0,60}},
                  lineColor={0,127,255},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}));
        end FluidNeutral;

      end Domains;

      package Junctions
        "<html><a href=\"modelica://Modelica\">Modelica</a> junctions between pure substances and their mixtures</html>"
        extends Modelica.Icons.Package;

        model Junction2
          "Junction between two pure substances and their mixture"
          import assert = FCSys.Utilities.assertEval;
          extends PartialJunction;

          replaceable package Medium1 =
              Modelica.Media.IdealGases.SingleGases.H2 constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium2 =
              Modelica.Media.IdealGases.SingleGases.H2O constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final
              package Medium = Medium1) "Fluid port for the 1st pure substance"
            annotation (Placement(transformation(extent={{30,30},{50,50}}),
                iconTransformation(extent={{30,30},{50,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final
              package Medium = Medium2) "Fluid port for the 2nd pure substance"
            annotation (Placement(transformation(extent={{30,-50},{50,-30}}),
                iconTransformation(extent={{30,-50},{50,-30}})));

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
            Documentation(info="<html><p>Assumptions:</p>
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol></html>"),
            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics={Line(
                  points={{0,-40},{40,-40}},
                  color={0,127,255},
                  smooth=Smooth.None), Line(
                  points={{0,40},{40,40}},
                  color={0,127,255},
                  smooth=Smooth.None)}));
        end Junction2;

        model Junction3
          "Junction between three pure substances and their mixture"
          import assert = FCSys.Utilities.assertEval;
          extends PartialJunction(redeclare replaceable package MixtureMedium
              = Media.CathodeGas);

          replaceable package Medium1 =
              Modelica.Media.IdealGases.SingleGases.H2O constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium2 =
              Modelica.Media.IdealGases.SingleGases.N2 constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));
          replaceable package Medium3 =
              Modelica.Media.IdealGases.SingleGases.O2 constrainedby
            Modelica.Media.Interfaces.PartialPureSubstance
            "<html>Medium model for the 3<sup>rd</sup> pure substance</html>"
            annotation (choicesAllMatching=true, Dialog(group=
                  "Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final
              package Medium = Medium1) "Fluid port for the 1st pure substance"
            annotation (Placement(transformation(extent={{30,30},{50,50}}),
                iconTransformation(extent={{30,30},{50,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final
              package Medium = Medium2) "Fluid port for the 2nd pure substance"
            annotation (Placement(transformation(extent={{30,-10},{50,10}}),
                iconTransformation(extent={{30,-10},{50,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b purePort3(redeclare final
              package Medium = Medium3) "Fluid port for the 3rd pure substance"
            annotation (Placement(transformation(extent={{30,-50},{50,-30}}),
                iconTransformation(extent={{30,-50},{50,-30}})));

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
            Documentation(info="<html><p>Assumptions:</p>
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol></html>"),
            Icon(graphics={
                Line(
                  points={{0,-40},{40,-40}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Line(
                  points={{0,40},{40,40}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Line(
                  points={{6,0},{40,0}},
                  color={0,127,255},
                  smooth=Smooth.None)}));
        end Junction3;

      protected
        partial model PartialJunction
          "Base model for a junction between pure substances and their mixture"
          extends FCSys.Icons.Names.Top3;

          replaceable package MixtureMedium = Media.AnodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium
            "Medium model for the mixture" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_a mixturePort(redeclare final
              package Medium = MixtureMedium) "Fluid port for the mixture"
            annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
                iconTransformation(extent={{-50,-10},{-30,10}})));

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
          annotation (defaultComponentName="junction", Icon(graphics={
                Line(
                  points={{-40,0},{0,0}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Line(
                  points={{0,-40},{0,40}},
                  color={0,127,255},
                  smooth=Smooth.None),
                Ellipse(
                  extent={{-6,6},{6,-6}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}));
        end PartialJunction;

      end Junctions;

      package Media
        "<html><a href=\"modelica://Modelica.Media\">Modelica.Media</a> models to interface with the <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">cell</a></html>"
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

    extends Modelica.Icons.Package;

    package BoundaryBus
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"
      extends Modelica.Icons.Package;

      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connectors</html>"

        extends Modelica.Icons.Package;

        model Temperature
          "<html>Impose a temperature difference across and a current between a pair of <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connectors</html>"

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

          Connectors.BoundaryBus negative
            "Negative-side multi-species connector for material, momentum, and energy"
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={-100,0})));
          Connectors.BoundaryBus positive
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
          annotation (Diagram(graphics));
        end Temperature;

        model Pressure
          "<html>Impose a pressure difference across and a heat flow rate between a pair of <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connectors</html>"
          extends Temperature(
            gas(
              H2(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature),
              H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature),
              N2(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature),
              O2(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature)),
            graphite(redeclare replaceable
                Conditions.ByConnector.ThermalDiffusive.Pair.HeatRate 'C+',
                'e-'(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature)),
            ionomer(
              redeclare replaceable
                Conditions.ByConnector.ThermalDiffusive.Pair.HeatRate 'SO3-',
              'H+'(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature),
              H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature)),
            liquid(H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Pair.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Pair.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Pair.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Pair.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Pair.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Pair.Thermal.temperature)));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          annotation (Diagram(graphics));

        end Pressure;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"
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

            Boundary.Pair.Temperature H2 if inclH2 "H2 conditions" annotation (
                Dialog(
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

            Boundary.Pair.Temperature H2O if inclH2O "H2O conditions"
              annotation (Dialog(
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

            Boundary.Pair.Temperature N2 if inclN2 "N2 conditions" annotation (
                Dialog(
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

            Boundary.Pair.Temperature O2 if inclO2 "O2 conditions" annotation (
                Dialog(
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

            ThermalDiffusive.Pair.HeatRate 'C+' if 'inclC+' "C+ conditions"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>C<sup>+</sup> conditions</html>",
                enable='inclC+'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

            parameter Boolean 'incle-'=false "Include e-" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Boundary.Pair.Temperature 'e-' if 'incle-' "e- conditions"
              annotation (Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>e<sup>-</sup> conditions</html>",
                enable='incle-'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

          equation
            // C+
            connect('C+'.negative, negative.'C+') annotation (Line(
                points={{-10,-0.2},{-10,5.55112e-016},{-100,5.55112e-016}},
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
                points={{6.10623e-016,-5},{-4.87687e-022,-50},{5.55112e-016,-50}},

                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // e-
            connect('e-'.negative, negative.'e-') annotation (Line(
                points={{-10,0},{-10,5.55112e-016},{-100,5.55112e-016}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('e-'.positive, positive.'e-') annotation (Line(
                points={{10,0},{100,0}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(u.'e-', 'e-'.u) annotation (Line(
                points={{5.55112e-016,50},{0,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('e-'.y, y.'e-') annotation (Line(
                points={{0,-5},{0,-50},{5.55112e-016,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                    extent={{-100,-100},{100,100}}), graphics));
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

            ThermalDiffusive.Pair.HeatRate 'SO3-' if 'inclSO3-'
              "SO3- conditions" annotation (Dialog(
                group="Species",
                __Dymola_label=
                    "<html>SO<sub>3</sub><sup>-</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='inclSO3-'), Placement(transformation(extent={{-10,-10},
                      {10,10}})));

            parameter Boolean 'inclH+'=false "Include H+" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            Boundary.Pair.Temperature 'H+' if 'inclH+' "H+ conditions"
              annotation (Dialog(
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

            Boundary.Pair.Temperature H2O if inclH2O "H2O conditions"
              annotation (Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // C19HF37O5S-
            connect('SO3-'.negative, negative.'SO3-') annotation (Line(
                points={{-10,-0.2},{-10,5.55112e-016},{-100,5.55112e-016}},
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
                points={{-10,0},{-10,5.55112e-016},{-100,5.55112e-016}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect('H+'.positive, positive.'H+') annotation (Line(
                points={{10,0},{10,5.55112e-016},{100,5.55112e-016}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(u.'H+', 'H+'.u) annotation (Line(
                points={{5.55112e-016,50},{0,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect('H+'.y, y.'H+') annotation (Line(
                points={{0,-5},{5.55112e-016,-50}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            // H2O
            connect(H2O.negative, negative.H2O) annotation (Line(
                points={{-10,0},{-10,5.55112e-016},{-100,5.55112e-016}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(H2O.positive, positive.H2O) annotation (Line(
                points={{10,0},{10,5.55112e-016},{100,5.55112e-016}},
                color={127,127,127},
                pattern=LinePattern.None,
                smooth=Smooth.None));

            connect(u.H2O, H2O.u) annotation (Line(
                points={{5.55112e-016,50},{0,5}},
                color={0,0,127},
                thickness=0.5,
                smooth=Smooth.None));

            connect(H2O.y, y.H2O) annotation (Line(
                points={{0,-5},{5.55112e-016,-50}},
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

            Boundary.Pair.Temperature H2O if inclH2O "H2O conditions"
              annotation (Dialog(
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

            Connectors.BoundaryBus negative
              "Negative-side multi-species connector for material, momentum, and energy"
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                  iconTransformation(extent={{-110,-10},{-90,10}})));
            Connectors.BoundaryBus positive
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
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"

        extends Modelica.Icons.Package;

        model Source
          "<html>Material source for a <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"

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

          Connectors.BoundaryBus boundary
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
          connect(gas.boundary, boundary.gas) annotation (Line(
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
          connect(graphite.boundary, boundary.graphite) annotation (Line(
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
          connect(ionomer.boundary, boundary.ionomer) annotation (Line(
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
          connect(liquid.boundary, boundary.liquid) annotation (Line(
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

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}), graphics));
        end Source;

        model Sink
          "<html>Material sink for a <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"
          extends Source(
            gas(
              H2(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0)),

              H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0)),

              N2(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0)),

              O2(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0))),

            graphite(redeclare replaceable ThermalDiffusive.Single.HeatRate
                'C+', 'e-'(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0))),

            ionomer(
              redeclare replaceable ThermalDiffusive.Single.HeatRate 'SO3-',
              'H+'(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0)),

              H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0))),

            liquid(H2O(
                redeclare replaceable function materialSpec =
                    Boundary.Single.Material.pressure,
                redeclare replaceable function afterSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function beforeSpec =
                    Boundary.Single.Translational.force,
                redeclare replaceable function thermalSpec =
                    Boundary.Single.Thermal.heatRate,
                redeclare replaceable function materialMeas =
                    Boundary.Single.Material.current,
                redeclare replaceable function afterMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function beforeMeas =
                    Boundary.Single.Translational.velocity,
                redeclare replaceable function thermalMeas =
                    Boundary.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=
                      U.atm),
                redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0))));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          // The materialSet and thermalSet blocks are redeclared as not replaceable
          // because y is set directly and can't be undone at instantiation.

          annotation (Diagram(graphics));

        end Sink;

        package Phases
          "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a> connector</html>"
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

            Boundary.Single.Source H2 if inclH2 "Include H2" annotation (Dialog(
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

            Boundary.Single.Source H2O if inclH2O "H2O conditions" annotation (
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

            Boundary.Single.Source N2 if inclN2 "N2 conditions" annotation (
                Dialog(
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

            Boundary.Single.Source O2 if inclO2 "O2 conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_label="<html>O<sub>2</sub> conditions</html>",
                enable=inclO2), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2
            connect(H2.boundary, boundary.H2) annotation (Line(
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
            connect(H2O.boundary, boundary.H2O) annotation (Line(
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
            connect(N2.boundary, boundary.N2) annotation (Line(
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
            connect(O2.boundary, boundary.O2) annotation (Line(
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

            replaceable
              Conditions.ByConnector.ThermalDiffusive.Single.Temperature 'C+'
              if 'inclC+' constrainedby
              Conditions.ByConnector.ThermalDiffusive.Single.Partial
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

            Boundary.Single.Source 'e-' if 'incle-' "e- conditions" annotation
              (Dialog(
                group="Species",
                __Dymola_label="<html>e<sup>-</sup> conditions</html>",
                __Dymola_descriptionLabel=true,
                enable='incle-'), Placement(transformation(extent={{-10,-10},{
                      10,10}})));

          equation
            // C+
            connect('C+'.therm, boundary.'C+') annotation (Line(
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
            connect('e-'.boundary, boundary.'e-') annotation (Line(
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
              "Include C19HF37O5S- (abbreviated as SO3-)" annotation (
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_label=
                    "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>, abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true));

            replaceable Conditions.ByConnector.ThermalDiffusive.Single.HeatRate
              'SO3-' if 'inclSO3-' constrainedby
              Conditions.ByConnector.ThermalDiffusive.Single.Partial
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

            Boundary.Single.Source 'H+' if 'inclH+' "H+ conditions" annotation
              (Dialog(
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

            Boundary.Single.Source H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // SO3-
            connect('SO3-'.therm, boundary.'SO3-') annotation (Line(
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
            connect('H+'.boundary, boundary.'H+') annotation (Line(
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
            connect(H2O.boundary, boundary.H2O) annotation (Line(
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

            Boundary.Single.Source H2O if inclH2O "H2O conditions" annotation (
                Dialog(
                group="Species",
                __Dymola_label="<html>H<sub>2</sub>O conditions</html>",
                __Dymola_descriptionLabel=true,
                enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,
                      10}})));

          equation
            // H2O
            connect(H2O.boundary, boundary.H2O) annotation (Line(
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

            Connectors.BoundaryBus boundary
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

    end BoundaryBus;

    package Boundary
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connector</html>"
      extends Modelica.Icons.Package;

      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connectors</html>"
        extends Modelica.Icons.Package;
        model Temperature
          "<html>Impose a temperature difference across and a current between a pair of <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connectors</html>"
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
          replaceable function afterSpec = Translational.velocity
            constrainedby Translational.Partial "Quantity" annotation (
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
          replaceable function beforeSpec = Translational.velocity
            constrainedby Translational.Partial "Quantity" annotation (
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
          replaceable function thermalSpec = Thermal.temperature constrainedby
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
          replaceable function afterMeas = Translational.force constrainedby
            Translational.Partial "First transverse quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // 2nd transverse
          replaceable function beforeMeas = Translational.force constrainedby
            Translational.Partial "Second transverse quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // Thermal
          replaceable function thermalMeas = Thermal.heatRate constrainedby
            Thermal.Partial "Thermal quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // Aliases
          Q.Pressure Deltap "Difference in pressure";
          Q.Velocity Deltaphi[Orient] "Difference in velocity";
          Q.Temperature DeltaT "Difference in temperature";

          Connectors.Boundary negative "Negative boundary" annotation (
              Placement(transformation(extent={{-110,-10},{-90,10}})));
          Connectors.Boundary positive "Positive boundary"
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

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}), graphics));
        end Temperature;

        model Pressure
          "<html>Impose a pressure difference across and a heat flow rate between a pair of <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connectors</html>"
          extends Temperature(
            redeclare replaceable function materialSpec = Material.pressure,
            redeclare replaceable function afterSpec = Translational.force,
            redeclare replaceable function beforeSpec = Translational.force,
            redeclare replaceable function thermalSpec = Thermal.heatRate,
            redeclare replaceable function materialMeas = Material.current,
            redeclare replaceable function afterMeas = Translational.velocity,
            redeclare replaceable function beforeMeas = Translational.velocity,

            redeclare replaceable function thermalMeas = Thermal.temperature);

          // See note in Reaction.Efforts.

        end Pressure;

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
            input Q.Pressure Deltap "Difference in pressure" annotation (Dialog(
                  __Dymola_label="<html>&Delta;<i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity Deltaphi[Orient] "Difference in velocity"
              annotation (Dialog(__Dymola_label="<html>&Delta;&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.Temperature DeltaT "Difference in temperature" annotation (
                Dialog(__Dymola_label="<html>&Delta;<i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
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
            input Q.Pressure Deltap "Difference in pressure" annotation (Dialog(
                  __Dymola_label="<html>&Delta;<i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity Deltaphi[Orient] "Difference in velocity"
              annotation (Dialog(__Dymola_label="<html>&Delta;&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.Temperature DeltaT "Difference in temperature" annotation (
                Dialog(__Dymola_label="<html>&Delta;<i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            input Orient orient
              "Orientation of translational momentum w.r.t. the boundary";

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));

            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

          end Partial;

        end Translational;

        package Thermal "Conditions for thermal diffusion"
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
            input Q.Pressure Deltap "Difference in pressure" annotation (Dialog(
                  __Dymola_label="<html>&Delta;<i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity Deltaphi[Orient] "Difference in velocity"
              annotation (Dialog(__Dymola_label="<html>&Delta;&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.Temperature DeltaT "Difference in temperature" annotation (
                Dialog(__Dymola_label="<html>&Delta;<i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

          end Partial;

        end Thermal;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connector</html>"
        extends Modelica.Icons.Package;
        model Source
          "<html>Material source for a <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connector</html>"
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
          replaceable function afterSpec = Translational.velocity
            constrainedby Translational.Partial "Quantity" annotation (
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
          replaceable function beforeSpec = Translational.velocity
            constrainedby Translational.Partial "Quantity" annotation (
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
          replaceable function thermalSpec = Thermal.temperature constrainedby
            Thermal.Partial "Quantity" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(tab="Specification", group="Thermal"),
            Placement(transformation(extent={{4,-10},{24,10}})));

          parameter Boolean internalThermal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(tab="Specification", group="Thermal"));
          replaceable Sources.RealExpression thermalSet(y=298.15*U.K) if
            internalThermal constrainedby Modelica.Blocks.Interfaces.SO
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
          // Aliases
          Q.PressureAbsolute p "Pressure";
          Q.TemperatureAbsolute T "Temperature";
          // Material
          replaceable function materialMeas = Material.pressure constrainedby
            Material.Partial "Material quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // 1st transverse
          replaceable function afterMeas = Translational.force constrainedby
            Translational.Partial "First transverse quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // 2nd transverse
          replaceable function beforeMeas = Translational.force constrainedby
            Translational.Partial "Second transverse quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          // Thermal
          replaceable function thermalMeas = Thermal.heatRate constrainedby
            Thermal.Partial "Thermal quantity" annotation (
              __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

          Connectors.Boundary boundary
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
                      boundary.p,
                      boundary.Ndot,
                      boundary.phi,
                      boundary.mPhidot,
                      boundary.T,
                      boundary.Qdot,
                      orient=Orient.after)
            "Internal, working value of first transverse specification"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={4,20})));
          Connectors.RealOutputInternal _u_before=beforeSpec(
                      boundary.p,
                      boundary.Ndot,
                      boundary.phi,
                      boundary.mPhidot,
                      boundary.T,
                      boundary.Qdot,
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
                        boundary.p,
                        boundary.Ndot,
                        boundary.phi,
                        boundary.mPhidot,
                        boundary.T,
                        boundary.Qdot)) "Generate the material output"
            annotation (Placement(transformation(extent={{40,50},{60,70}})));
          Sources.RealExpression beforeOut(y=beforeMeas(
                        boundary.p,
                        boundary.Ndot,
                        boundary.phi,
                        boundary.mPhidot,
                        boundary.T,
                        boundary.Qdot,
                        orient=Orient.before))
            "Generate the 2nd transverse output"
            annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
          Sources.RealExpression afterOut(y=afterMeas(
                        boundary.p,
                        boundary.Ndot,
                        boundary.phi,
                        boundary.mPhidot,
                        boundary.T,
                        boundary.Qdot,
                        orient=Orient.after))
            "Generate the 1st transverse output"
            annotation (Placement(transformation(extent={{40,10},{60,30}})));
          Sources.RealExpression thermalOut(y=thermalMeas(
                        boundary.p,
                        boundary.Ndot,
                        boundary.phi,
                        boundary.mPhidot,
                        boundary.T,
                        boundary.Qdot)) "Generate the thermal output"
            annotation (Placement(transformation(extent={{40,-70},{60,-50}})));
        equation
          // Aliases
          p = boundary.p;
          T = boundary.T;

          _u_material = materialSpec(
                    boundary.p,
                    boundary.Ndot,
                    boundary.phi,
                    boundary.mPhidot,
                    boundary.T,
                    boundary.Qdot);
          _u_thermal = thermalSpec(
                    boundary.p,
                    boundary.Ndot,
                    boundary.phi,
                    boundary.mPhidot,
                    boundary.T,
                    boundary.Qdot);

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
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}), graphics));
        end Source;

        model Sink
          "<html>Material sink for a <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connector</html>"
          extends Source(
            redeclare replaceable function materialSpec = Material.pressure,
            redeclare replaceable function afterSpec = Translational.force,
            redeclare replaceable function beforeSpec = Translational.force,
            redeclare replaceable function thermalSpec = Thermal.heatRate,
            redeclare replaceable function materialMeas = Material.current,
            redeclare replaceable function afterMeas = Translational.velocity,
            redeclare replaceable function beforeMeas = Translational.velocity,

            redeclare replaceable function thermalMeas = Thermal.temperature,
            redeclare Modelica.Blocks.Sources.RealExpression materialSet(y=U.atm),

            redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=0));

          // Note:  In Dymola 2014, the paths must be explicitly given to prevent
          // the error "Cannot show the parameter dialog for redeclared class [...]".

          // The materialSet and thermalSet blocks are redeclared as not replaceable
          // because y is set directly and can't be undone at instantiation.

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}), graphics));

        end Sink;

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

          function volumeRate "Volumetric flow rate"
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

          function standardVolumeRate "Standard volumetric flow rate"
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
            parameter Q.TemperatureAbsolute To=298.15*U.K
              "Standard temperature" annotation (Dialog(__Dymola_label=
                    "<html><i>T</i><sup>o</sup></html>"));
            parameter Q.PressureAbsolute pp=U.bar "Standard pressure"
              annotation (Dialog(__Dymola_label=
                    "<html><i>p</i><sup>o</sup></html>"));

          algorithm
            x := Data.v_Tp(To, p0)*Ndot;
            annotation (Inline=true);
          end standardVolumeRate;

          partial function Partial
            "Template of a function to select a material quantity"
            extends Modelica.Icons.Function;

            // Material
            input Q.PressureAbsolute p "Pressure"
              annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity phi[Orient] "Velocity"
              annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.TemperatureAbsolute T "Temperature"
              annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

          end Partial;

          function potential "Gibbs potential"
            extends Partial;

            replaceable package Data =
                Characteristics.BaseClasses.Characteristic constrainedby
              Characteristics.BaseClasses.Characteristic "Characteristic data"
              annotation (
              Dialog(group="Material properties"),
              choicesAllMatching=true,
              __Dymola_choicesFromPackage=true,
              Placement(transformation(extent={{-60,40},{-40,60}}),
                  iconTransformation(extent={{-10,90},{10,110}})));

          algorithm
            x := Data.g(T, p);
            annotation (Inline=true);
          end potential;

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
            input Q.PressureAbsolute p "Pressure"
              annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity phi[Orient] "Velocity"
              annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.TemperatureAbsolute T "Temperature"
              annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            input Orient orient
              "Orientation of translational momentum w.r.t. the boundary";

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));

            annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

          end Partial;

        end Translational;

        package Thermal "Conditions for thermal diffusion"
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
            input Q.PressureAbsolute p "Pressure"
              annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
            input Q.Current Ndot "Current"
              annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));

            // Translational
            input Q.Velocity phi[Orient] "Velocity"
              annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
            input Q.Force mPhidot[Orient] "Shear force" annotation (Dialog(
                  __Dymola_label="<html><i>m</i>&Phi;dot</html>"));

            // Thermal
            input Q.TemperatureAbsolute T "Temperature"
              annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
            input Q.Power Qdot "Rate of thermal conduction"
              annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

            output Real x "Value of condition"
              annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
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

    end Boundary;

    package Amagat
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector</html>"
      extends Modelica.Icons.Package;

      model Pressure "Specify pressure (measure volume)"
        extends FCSys.Conditions.ByConnector.Amagat.Partial(final y=amagat.V,
            set(y=U.atm));

      equation
        amagat.p = u_final;

      end Pressure;

      model Volume "Provide volume (measure pressure)"
        extends Partial(final y=amagat.p, set(y=U.cc));

      equation
        amagat.V = u_final;

      end Volume;

      model VolumeFixed "Model to establish a fixed total volume"
        extends FCSys.Icons.Names.Top3;

        parameter Q.Volume V "Volume"
          annotation (Dialog(__Dymola_label="<html><i>V</i></html>"));

        Connectors.Amagat amagat "Connector for additivity of volume"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
              iconTransformation(extent={{-10,-10},{10,10}})));

      equation
        amagat.V = V;

        annotation (
          defaultComponentName="volume",
          Documentation(info="<html><p>This model uses an <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector that imposes
    additivity of volume.  In order to use additivity of pressure, use
    the <a href=\"modelica://FCSys.Conditions.Adapters.AmagatDalton\">AmagatDalton</a> adapter.</p>

    <p>See also the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

          Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
                  100,100}}), graphics={Polygon(
                      points={{-60,-60},{-60,20},{-20,60},{60,60},{60,-20},{20,
                  -60},{-60,-60}},
                      lineColor={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      fillColor={225,225,225},
                      fillPattern=FillPattern.Solid),Line(
                      points={{60,60},{20,20}},
                      color={0,0,0},
                      pattern=LinePattern.Dash,
                      smooth=Smooth.None),Line(
                      points={{-60,20},{20,20},{20,-60}},
                      color={0,0,0},
                      pattern=LinePattern.Dash,
                      smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                  {100,100}}), graphics));
      end VolumeFixed;

      partial model Partial "Base model for a pressure/volume"

        extends FCSys.Icons.Conditions.SingleShort;

        parameter Boolean internal=true "Use internal specification"
          annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Specification"));

        replaceable Modelica.Blocks.Sources.RealExpression set if internal
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

        connect(set.y, u_final) annotation (Line(
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
                  fillColor={255,255,255}),Ellipse(
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
        extends Partial(final y=dalton.p, set(y=U.cc));

      equation
        dalton.V = u_final;

      end Volume;

      model Pressure "Specify pressure (measure volume)"
        extends Partial(final y=dalton.V, set(y=U.atm));

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

        replaceable Modelica.Blocks.Sources.RealExpression set if internal
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

        connect(set.y, u_final) annotation (Line(
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

    package Inter
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> or <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a> connector</html>"
      extends Modelica.Icons.Package;

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
        replaceable function thermalSpec = ThermalDiffusive.heatRate
          constrainedby ThermalDiffusive.Partial "Quantity" annotation (
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
            Conditions.ByConnector.Inter.Translational.velocity constrainedby
          Translational.Partial "X-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));
        //
        // Y-axis translational
        replaceable function transYMeas =
            Conditions.ByConnector.Inter.Translational.velocity constrainedby
          Translational.Partial "Y-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));
        //
        // Z-axis translational
        replaceable function transZMeas =
            Conditions.ByConnector.Inter.Translational.velocity constrainedby
          Translational.Partial "Z-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));
        //
        // Thermal
        replaceable function thermalMeas =
            Conditions.ByConnector.Inter.ThermalDiffusive.temperature
          constrainedby ThermalDiffusive.Partial "Thermal quantity" annotation
          (__Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

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

        Connectors.Inert inter(final n_trans=n_trans)
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
          "X-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        final Connectors.RealOutput y_transY=transYMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,20}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));

        final Connectors.RealOutput y_transZ=transZMeas(
                  inter.phi,
                  inter.mPhidot,
                  inter.T,
                  inter.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Placement(
              transformation(
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
                  inter.Qdot) "Thermal measurement" annotation (Placement(
              transformation(
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

      model Efforts
        "<html>Condition for an <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a> or <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> connector, with efforts specified by default</html>"

        extends FCSys.Conditions.ByConnector.Inter.Flows(
          redeclare replaceable function transXSpec = Translational.velocity,
          redeclare replaceable function transYSpec = Translational.velocity,
          redeclare replaceable function transZSpec = Translational.velocity,
          redeclare replaceable function thermalSpec =
              ThermalDiffusive.temperature,
          redeclare replaceable function transXMeas = Translational.force,
          redeclare replaceable function transYMeas = Translational.force,
          redeclare replaceable function transZMeas = Translational.force,
          redeclare replaceable function thermalMeas =
              ThermalDiffusive.heatRate,
          redeclare Modelica.Blocks.Sources.RealExpression thermalSet(y=298.15*
                U.K));

        // The thermalSet block is redeclared as not replaceable because
        // y is set directly and can't be undone at instantiation.

        annotation (defaultComponentName="inert", Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics));

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
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          // Thermal
          input Q.TemperatureAbsolute T "Temperature"
            annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
          input Q.Power Qdot "Rate of thermal conduction"
            annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

          input Integer i(min=1,max=3) "Index of the translational axis"
            annotation (Dialog(__Dymola_label="<html><i>i</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

        end Partial;

      end Translational;

      package ThermalDiffusive "Conditions for thermal diffusion"
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
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          // Thermal
          input Q.TemperatureAbsolute T "Temperature"
            annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
          input Q.Power Qdot "Rate of thermal conduction"
            annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

        end Partial;

      end ThermalDiffusive;
      annotation (Icon(graphics={Ellipse(
                  extent={{-60,60},{60,-60}},
                  lineColor={170,0,0},
                  fillPattern=FillPattern.Solid,
                  fillColor={221,23,47})}));

    end Inter;

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
                  trans.phi,
                  trans.mPhidot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));

        final Connectors.RealOutput y_transY=transYMeas(
                  trans.phi,
                  trans.mPhidot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));

        final Connectors.RealOutput y_transZ=transZMeas(
                  trans.phi,
                  trans.mPhidot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40})));

        Connectors.Translational trans(final n_trans=n_trans)
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
                  trans.phi,
                  trans.mPhidot,
                  i=transCart[Axis.x]) if inclTransX
          "Internal, working value of X-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,40})));

        Connectors.RealOutputInternal _u_transY=transYSpec(
                  trans.phi,
                  trans.mPhidot,
                  i=transCart[Axis.y]) if inclTransY
          "Internal, working value of Y-axis translational specification"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-36,0})));

        Connectors.RealOutputInternal _u_transZ=transZSpec(
                  trans.phi,
                  trans.mPhidot,
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

        annotation (Icon(graphics));
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
          // Translational
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          input Integer i(min=1,max=3) "Index of the translational axis"
            annotation (Dialog(__Dymola_label="<html><i>i</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
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

    package ThermalDiffusive
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.ThermalDiffusive\">ThermalDiffusive</a> connector</html>"
      extends Modelica.Icons.Package;

      package Pair
        "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.ThermalDiffusive\">ThermalDiffusive</a> connectors</html>"
        extends Modelica.Icons.Package;

        model HeatRate "Specify heat flow rate (measure temperature)"
          extends Partial(final y=positive.T - negative.T);

        equation
          negative.Qdot = u_final;

        end HeatRate;

        model Temperature "Specify temperature (measure heat flow rate)"
          extends Partial(final y=negative.Qdot);

        equation
          positive.T - negative.T = u_final;

        end Temperature;

        partial model Partial "Base model for a thermal condition"

          extends FCSys.Icons.Conditions.PairShort;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression set if internal
            constrainedby Modelica.Blocks.Interfaces.SO
            "Source of internal specification" annotation (
            __Dymola_choicesFromPackage=true,
            Dialog(group="Specification",enable=internal),
            Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={10,40})));

          Connectors.RealInput u if not internal "Value of specified condition"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,110}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,50})));

          Connectors.RealOutput y "Measurement expression" annotation (Dialog(
                tab="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-110}),iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-50})));

          Connectors.ThermalDiffusive negative
            "Negative connector for thermal diffusion"
            annotation (Placement(transformation(extent={{-110,-12},{-90,8}})));

        protected
          Connectors.RealOutputInternal u_final
            "Final value of specified condition" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,8})));

        public
          Connectors.ThermalDiffusive positive
            "Positive connector for thermal diffusion"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));
        equation

          // Conservation of energy
          0 = negative.Qdot + positive.Qdot;
          connect(u, u_final) annotation (Line(
              points={{0,110},{0,8}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(set.y, u_final) annotation (Line(
              points={{10,29},{10,20},{0,20},{0,8}},
              color={0,0,127},
              smooth=Smooth.None));

          annotation (
            defaultComponentName="thermal",
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics),
            Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                    {100,100}}), graphics));
        end Partial;

      end Pair;

      package Single
        "<html>Conditions for a single <a href=\"modelica://FCSys.Connectors.ThermalDiffusive\">ThermalDiffusive</a> connector</html>"
        extends Modelica.Icons.Package;

        model HeatRate "Specify heat flow rate (measure temperature)"
          extends Partial(final y=therm.T);

        equation
          therm.Qdot = u_final;

        end HeatRate;

        model Temperature "Specify temperature (measure heat flow rate)"
          extends Partial(final y=therm.Qdot, set(y=298.15*U.K));

        equation
          therm.T = u_final;

        end Temperature;

        model Resistance
          "Specify thermal resistance to the environment (measure temperature)"
          extends Partial(final y=therm.T);

          outer Conditions.Environment environment "Environmental conditions";

        equation
          therm.T - environment.T = u_final*therm.Qdot;

        end Resistance;

        partial model Partial "Base model for a thermal condition"

          extends FCSys.Icons.Conditions.SingleShort;

          parameter Boolean internal=true "Use internal specification"
            annotation (
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(group="Specification"));

          replaceable Modelica.Blocks.Sources.RealExpression set if internal
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
                tab="Measurement"), Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={110,0})));

          Connectors.ThermalDiffusive therm "Connector for thermal diffusion"
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

          connect(set.y, u_final) annotation (Line(
              points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
              color={0,0,127},
              smooth=Smooth.None));

          annotation (defaultComponentName="thermal");
        end Partial;

      end Single;
      annotation (Icon(graphics={Ellipse(
                  extent={{-60,60},{60,-60}},
                  lineColor={127,127,127},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,255,255})}));

    end ThermalDiffusive;

    package Chemical
      "<html>Conditions for an <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connector</html>"
      extends Modelica.Icons.Package;

      model Potential "Specify chemical potential (measure current)"
        extends Partial(final y=chemical.Ndot);

      equation
        chemical.g = u_final;

      end Potential;

      model Current "Specify current (measure chemical potential)"
        extends Partial(final y=chemical.g);

      equation
        chemical.Ndot = u_final;

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

        replaceable Modelica.Blocks.Sources.RealExpression set if internal
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
          "Product of specific entropy and temperature" annotation (Dialog(
              group="Properties upon outflow", __Dymola_label=
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
        output Q.Velocity phi_actual[n_trans]=actualStream(chemical.phi)
          "Velocity of the actual stream";
        output Q.Potential sT_actual=actualStream(chemical.sT)
          "Specific entropy-temperature product of the actual stream";
        Connectors.Chemical chemical(n_trans=n_trans)
          "Connector for a species of a chemical reaction"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

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
        chemical.phi = phi[cartTrans];
        chemical.sT = sT;

        connect(set.y, u_final) annotation (Line(
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
                  lineColor={239,142,1},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,195,38})}));

    end Chemical;

    package Reaction
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector</html>"
      extends Modelica.Icons.Package;

      model Rate
        "<html>Impose a reaction rate via a <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector</html>"
        import Modelica.Math.BooleanVectors.countTrue;
        import Modelica.Math.BooleanVectors.enumerate;
        import Modelica.Blocks.Sources;
        extends FCSys.Icons.Conditions.Single;

        // Specification
        // -------------
        // Material
        replaceable function materialSpec = Material.reactionRate
          constrainedby Material.Partial "Quantity" annotation (
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
        replaceable function transXSpec = Translational.velocity constrainedby
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
        replaceable function transYSpec = Translational.velocity constrainedby
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
        replaceable function transZSpec = Translational.velocity constrainedby
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
        replaceable function thermalSpec =
            ThermalAdvective.specificEntropyTemperature constrainedby
          ThermalAdvective.Partial "Quantity" annotation (
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
        replaceable function materialMeas = Material.potential constrainedby
          Material.Partial "Material quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // X-axis translational
        replaceable function transXMeas = Translational.force constrainedby
          Translational.Partial "X-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Y-axis translational
        replaceable function transYMeas = Translational.force constrainedby
          Translational.Partial "Y-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Z-axis translational
        replaceable function transZMeas = Translational.force constrainedby
          Translational.Partial "Z-axis translational quantity" annotation (
            __Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

        // Thermal
        replaceable function thermalMeas = ThermalAdvective.heatRate
          constrainedby ThermalAdvective.Partial "Thermal quantity" annotation
          (__Dymola_choicesFromPackage=true, Dialog(tab="Measurement"));

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
                  reaction.g,
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
                  reaction.g,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.x]) if inclTransX
          "X-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0})));
        final Connectors.RealOutput y_transY=transYMeas(
                  reaction.Ndot,
                  reaction.g,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.y]) if inclTransY
          "Y-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,40})));
        final Connectors.RealOutput y_transZ=transZMeas(
                  reaction.Ndot,
                  reaction.g,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot,
                  i=transCart[Axis.z]) if inclTransZ
          "Z-axis translational measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-40})));
        final Connectors.RealOutput y_thermal=thermalMeas(
                  reaction.Ndot,
                  reaction.g,
                  reaction.phi,
                  reaction.mPhidot,
                  reaction.sT,
                  reaction.Qdot) "Thermal measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-80}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={110,-80})));

        Connectors.Reaction reaction(final n_trans=n_trans)
          "Stoichiometric connector for the chemical reaction"
          annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

      protected
        final inner parameter Integer n_trans=countTrue({inclTransX,inclTransY,
            inclTransZ}) "Number of components of translational momentum";
        final inner parameter Integer transCart[Axis]=enumerate({inclTransX,
            inclTransY,inclTransZ})
          "Translational-momentum-component indices of the Cartesian axes";

        Connectors.RealOutputInternal _u_material=materialSpec(
                  reaction.Ndot,
                  reaction.g,
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
                  reaction.g,
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
                  reaction.g,
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
                  reaction.g,
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
                  reaction.g,
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
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent=
                  {{-100,-100},{100,100}}), graphics));
      end Rate;

      model Offset
        "<html>Add a potential to a <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector</html>"

        extends FCSys.Conditions.ByConnector.Reaction.Rate(
          redeclare replaceable function materialSpec = Material.potential,
          redeclare replaceable function transXSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function transYSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function transZSpec =
              Conditions.ByConnector.Reaction.Translational.velocity,
          redeclare replaceable function thermalSpec =
              Conditions.ByConnector.Reaction.ThermalAdvective.specificEntropyTemperature,

          redeclare replaceable function materialMeas =
              Conditions.ByConnector.Reaction.Material.reactionRate,
          redeclare replaceable function transXMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function transYMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function transZMeas =
              Conditions.ByConnector.Reaction.Translational.force,
          redeclare replaceable function thermalMeas =
              Conditions.ByConnector.Reaction.ThermalAdvective.heatRate);

        // Note:  In Dymola 2014, the paths must be explicitly given to prevent
        // the error "Cannot show the parameter dialog for redeclared class [...]".

      end Offset;

      package Material "Material conditions"
        extends Modelica.Icons.Package;

        function reactionRate "Rate of the reaction"
          extends Partial;

        algorithm
          x := Ndot;
          annotation (Inline=true);
        end reactionRate;

        function potential "Offset from the potential of the reaction"
          extends Partial;

        algorithm
          x := g;
          annotation (Inline=true);
        end potential;

        partial function Partial
          "Template of a function to select a material quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "Rate of reaction"
            annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));
          input Q.Potential g "Chemical potential"
            annotation (Dialog(__Dymola_label="<html><i>g</i></html>"));

          // Translational advection
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          // Thermal advection
          input Q.PotentialAbsolute sT
            "Product of specific entropy and temperature"
            annotation (Dialog(__Dymola_label="<html><i>sT</i></html>"));
          input Q.Power Qdot "Rate of thermal advection"
            annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

        end Partial;

      end Material;

      package Translational "Translational conditions"
        extends Modelica.Icons.Package;

        function velocity "Velocity of the product"
          extends Partial;

        algorithm
          x := phi[i];
          annotation (Inline=true);
        end velocity;

        function force "Force from the stream"
          extends Partial;

        algorithm
          x := mPhidot[i];
          annotation (Inline=true);
        end force;

        partial function Partial
          "Template of a function to select a translational quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "Rate of reaction"
            annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));
          input Q.Potential g "Chemical potential"
            annotation (Dialog(__Dymola_label="<html><i>g</i></html>"));

          // Translational advection
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          // Thermal advection
          input Q.PotentialAbsolute sT
            "Product of specific entropy and temperature"
            annotation (Dialog(__Dymola_label="<html><i>sT</i></html>"));
          input Q.Power Qdot "Rate of thermal advection"
            annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

          input Integer i(min=1,max=3) "Index of the translational axis"
            annotation (Dialog(__Dymola_label="<html><i>i</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

        end Partial;

      end Translational;

      package ThermalAdvective "Conditions for thermal advection"
        extends Modelica.Icons.Package;

        function specificEntropyTemperature
          "Product of specific entropy and temperature of the product"
          extends Partial;

        algorithm
          x := sT;
          annotation (Inline=true);
        end specificEntropyTemperature;

        function heatRate "Rate of heat from the stream"
          extends Partial;

        algorithm
          x := Qdot;
          annotation (Inline=true);
        end heatRate;

        partial function Partial
          "Template of a function to select a thermal quantity"
          extends Modelica.Icons.Function;

          // Material diffusion
          input Q.Current Ndot "Rate of reaction"
            annotation (Dialog(__Dymola_label="<html><i>N&#775;</i></html>"));
          input Q.Potential g "Chemical potential"
            annotation (Dialog(__Dymola_label="<html><i>g</i></html>"));

          // Translational advection
          input Q.Velocity phi[:] "Velocity"
            annotation (Dialog(__Dymola_label="<html>&phi;</html>"));
          input Q.Force mPhidot[:] "Force"
            annotation (Dialog(__Dymola_label="<html><i>m</i>&Phi;dot</html>"));

          // Thermal advection
          input Q.PotentialAbsolute sT
            "Product of specific entropy and temperature"
            annotation (Dialog(__Dymola_label="<html><i>sT</i></html>"));
          input Q.Power Qdot "Rate of thermal advection"
            annotation (Dialog(__Dymola_label="<html><i>Q&#775;</i></html>"));

          output Real x "Value of condition"
            annotation (Dialog(__Dymola_label="<html><i>x</i></html>"));
          annotation (Inline=true, Documentation(info="<html>
  <p>This function takes as inputs all the efforts and flows of the associated
  connector.  It should be extended to add an algorithm that maps these inputs
  to a single value.</p></html>"));

        end Partial;

      end ThermalAdvective;

      annotation (Icon(graphics={Ellipse(
                  extent={{-60,60},{60,-60}},
                  lineColor={239,142,1},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,255,255}),Ellipse(
                  extent={{-30,30},{30,-30}},
                  fillColor={255,195,38},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0})}));

    end Reaction;
    annotation (Documentation(info="<html>
  <p>This package contains models to impose conditions on each of the declarative connectors
  established in <a href=\"modelica://FCSys.Connectors\">FCSys.Connectors</a>.  The subpackages
  are named according to the corresponding connector.</p>
</html>"));

  end ByConnector;

  record Environment "Environmental properties for a simulation"
    extends FCSys.Icons.Names.Top3;

    // Store the values of the base constants and units.
    final constant U.Bases.Base baseUnits=U.base "Base constants and units";

    parameter Boolean analysis=true "Include optional variables for analysis"
      annotation (choices(__Dymola_checkBox=true));

    // Thermodynamics
    parameter Q.TemperatureAbsolute T(nominal=300*U.K) = 298.15*U.K
      "Temperature" annotation (Dialog(__Dymola_label="<html><i>T</i></html>",
          group="Thermodynamics"));
    parameter Q.PressureAbsolute p(nominal=U.atm) = U.atm "Pressure"
      annotation (Dialog(__Dymola_label="<html><i>p</i></html>", group=
            "Thermodynamics"));
    parameter Q.NumberAbsolute RH(
      displayUnit="%",
      max=1) = 0.5 "Relative humidity"
      annotation (Dialog(group="Thermodynamics"));
    parameter Q.NumberAbsolute psi_O2_dry(
      final max=1,
      displayUnit="%") = 0.20946
      "<html>Mole fraction of O<sub>2</sub> in the dry gas</html>" annotation (
        Dialog(__Dymola_label="<html>&psi;<sub>O2 dry</sub></html>", group=
            "Thermodynamics"));
    // Value from http://en.wikipedia.org/wiki/Oxygen, accessed 2013/10/30
    final parameter Q.PressureAbsolute p_sat=Characteristics.H2O.p_sat(T)
      "Saturation pressure of H2O vapor";
    final parameter Q.PressureAbsolute p_H2O=RH*p_sat "Pressure of H2O vapor";
    final parameter Q.PressureAbsolute p_dry=p - p_H2O "Pressure of dry gases";
    final parameter Q.PressureAbsolute p_O2=psi_O2_dry*p_dry "Pressure of O2";
    final parameter Q.NumberAbsolute psi_H2O=p_H2O/p "Mole fraction of H2O";
    final parameter Q.NumberAbsolute psi_dry=1 - psi_H2O
      "Mole fraction of dry gases";

    // Fields
    parameter Q.Acceleration a[Axis]={0,Modelica.Constants.g_n*U.m/U.s^2,0}
      "Acceleration due to body forces" annotation (Dialog(__Dymola_label=
            "<html><b><i>a</i></b></html>", group="Fields"));
    // The gravity component is positive because it's added to the transient
    // term in the Species model.
    parameter Q.ForceSpecific E[Axis]={0,0,0} "Electric field" annotation (
        Dialog(__Dymola_label="<html><b><i>E</i></b></html>", group="Fields"));

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

  model Router "Connect two pairs of boundaries to pass through or cross over"
    extends FCSys.Icons.Names.Top3;
    parameter Boolean crossOver=false "Cross over (otherwise, pass through)"
      annotation (choices(__Dymola_checkBox=true));
    Connectors.BoundaryBus negative1 "Negative boundary 1" annotation (
        Placement(transformation(extent={{-90,-50},{-70,-30}}, rotation=0),
          iconTransformation(extent={{-90,-50},{-70,-30}})));
    Connectors.BoundaryBus positive1 "Positive boundary 1" annotation (
        Placement(transformation(extent={{70,-50},{90,-30}}, rotation=0),
          iconTransformation(extent={{70,-50},{90,-30}})));
    Connectors.BoundaryBus negative2 "Negative boundary 2" annotation (
        Placement(transformation(extent={{-90,30},{-70,50}}, rotation=0),
          iconTransformation(extent={{-90,30},{-70,50}})));
    Connectors.BoundaryBus positive2 "Positive boundary 2" annotation (
        Placement(transformation(extent={{70,30},{90,50}}, rotation=0),
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
            color={225,225,225},
            thickness=0.5,
            smooth=Smooth.Bezier,
            pattern=LinePattern.Dash),
          Line(
            points={{-80,-40},{-40,-40},{0,0},{40,40},{80,40}},
            color={225,225,225},
            thickness=0.5,
            smooth=Smooth.Bezier,
            pattern=LinePattern.Dash),
          Line(
            points={{-82,40},{78,40}},
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
            points={{-82,40},{78,40}},
            color={225,225,225},
            visible=crossOver,
            smooth=Smooth.None,
            pattern=LinePattern.Dash,
            thickness=0.5),
          Line(
            points={{-80,-40},{80,-40}},
            color={225,225,225},
            visible=crossOver,
            smooth=Smooth.None,
            thickness=0.5,
            pattern=LinePattern.Dash),
          Line(
            points={{-80,40},{-40,40},{0,0},{40,-40},{80,-40}},
            color={127,127,127},
            thickness=0.5,
            visible=crossOver,
            smooth=Smooth.Bezier),
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
