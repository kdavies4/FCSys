within FCSys;
package BCs "Models for boundary conditions"
  extends Modelica.Icons.SourcesPackage;

  // TODO:  Recheck this package, fix errors and warnings.

  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;
    model Environment "<html>Test the <code>Environment</code> model</html>"
      extends Modelica.Icons.Example;

      // TODO:  Make this into a meaningful example.
      FCSys.BCs.Defaults default
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    end Environment;

    model FaceBC "<html>Test the BCs for the face of a subregion</html>"
      extends Modelica.Icons.Example;

      Face.Subregion subregionFaceBC(gas(inclH2O=true, H2O(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic,
              materialSpec(k=-0.4805*U.V))))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));
      Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclReact=false,
        inclXFaces=false,
        inclYFaces=true,
        inclZFaces=false,
        inclVelX=false,
        inclVelY=true,
        graphite(inclC=true, C(V_IC=0.5*U.cm^3,alpha_tau=1e-3*U.cm/U.A)),
        gas(
          inclH2O=true,
          H2O(
            alpha_tau=1e-3*U.cm/U.A,
            xNegative(
              thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic,

              slipY=false,
              slipZ=false),
            xPositive(
              thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic,

              slipY=false,
              slipZ=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false),
            yPositive(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic),

            initMethPartNum=FCSys.Subregions.BaseClasses.InitMethScalar.PotentialElectrochemical,

            mu_IC=-298685),
          inclH2=false,
          inclO2=false))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Defaults defaults
        annotation (Placement(transformation(extent={{30,30},{50,50}})));
    equation
      connect(subregion.yPositive, subregionFaceBC.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=0.0015, NumberOfIntervals=5000),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.FaceBC.mos"));
    end FaceBC;

    model FaceBCPhases
      "<html>Test the BCs for the face of a subregion with phases</html>"
      extends Modelica.Icons.Example;

      FCSys.BCs.Face.Phases.Phase phaseFaceBC(inclH2O=true, H2O(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic,
            redeclare FCSys.BCs.Face.Species.Material.Current materialBC))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));
      FCSys.Subregions.Phases.Phase subregion(
        inclReact=false,
        inclH2=false,
        inclH2O=true,
        H2O(
          xNegative(
            thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic,
            slipY=false,
            slipZ=false),
          xPositive(
            thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic,
            slipY=false,
            slipZ=false),
          yPositive(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic),

          zNegative(slipX=false, slipY=false),
          zPositive(slipX=false, slipY=false)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Defaults defaults
        annotation (Placement(transformation(extent={{30,30},{50,50}})));
    equation
      connect(subregion.yPositive, phaseFaceBC.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=0.003),
        experimentSetupOutput);
    end FaceBCPhases;

    model Router "<html>Test the <code>Router<code> model</html>"
      extends Modelica.Icons.Example;

      // TODO:  Make this into a meaningful example.
      FCSys.BCs.Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    end Router;

    model Adapteminus "<html>Test the <code>'e-Adapt'</code> model</html>"

      extends Modelica.Icons.Example;

      FCSys.BCs.Adapters.'AdaptSubregione-' 'e-Adapt'
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{30,0},{50,20}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T
          =298.15)
        annotation (Placement(transformation(extent={{50,-30},{30,-10}})));
      Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclReact=false,
        inclYFaces=false,
        inclZFaces=false,
        gas(inclH2=false, inclH2O=false),
        graphite('incle-'=true,'e-'(
            xNegative(
              slipY=false,
              slipZ=false,
              thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic),

            xPositive(
              slipY=false,
              slipZ=false,
              thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic),
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false),
            mu_IC=0.01*U.V)))
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

      inner Defaults defaults(T=350*U.K)
        annotation (Placement(transformation(extent={{60,40},{80,60}})));
    equation
      connect(ground.p, 'e-Adapt'.pin) annotation (Line(
          points={{40,20},{20,20},{20,4},{10,4}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(fixedTemperature.port, 'e-Adapt'.port) annotation (Line(
          points={{30,-20},{20,-20},{20,-4},{10,-4}},
          color={191,0,0},
          smooth=Smooth.None));

      connect(subregion.xPositive, 'e-Adapt'.face) annotation (Line(
          points={{-30,6.10623e-16},{-20,-3.36456e-22},{-20,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=2e-10),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.'e-Adapt'.mos"));
    end Adapteminus;

    model AdaptFluid "<html>Test the <code>FluidAdapt</code> model</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      FCSys.BCs.Adapters.AdaptBusH2 H2Adapt(redeclare package Medium =
            Modelica.Media.IdealGases.SingleGases.H2)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Modelica.Fluid.Vessels.ClosedVolume volume(
        use_portsData=false,
        nPorts=1,
        redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2 (
              referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false),
        V=1e-6) annotation (Placement(transformation(extent={{20,10},{40,30}})));
      inner Modelica.Fluid.System system
        annotation (Placement(transformation(extent={{40,70},{60,90}})));
      FCSys.Subregions.FlatSubregion subregion(
        L={1,1,1}*U.cm,
        inclH2=true,
        inclH2O=false,
        inclReact=false,
        inclYFaces=false,
        inclZFaces=false,
        H2(
          xPositive(
            thermoOpt=ThermoOpt.OpenDiabatic,
            slipY=false,
            slipZ=false),
          initMethPartNum=FCSys.Subregions.BaseClasses.InitMethScalar.None,
          xNegative(
            final thermoOpt=ThermoOpt.ClosedAdiabatic,
            slipY=false,
            slipZ=false)))
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      inner Defaults defaults(analysis=true, T=293.15*U.K)
        annotation (Placement(transformation(extent={{70,70},{90,90}})));

    equation
      connect(H2Adapt.fluidPort, volume.ports[1]) annotation (Line(
          points={{10,6.10623e-16},{30,6.10623e-16},{30,10}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(subregion.xPositive, H2Adapt.face) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=1.5e-06),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.FluidAdapt.mos"));
    end AdaptFluid;

    model AdaptFluid2 "<html>Test the <code>FluidAdapt</code> model</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      ClosedVolume volume(
        use_portsData=false,
        redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2,
        V=1e-6)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
      inner Modelica.Fluid.System system
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
      FCSys.Subregions.FlatSubregion subregion(
        L={1,1,1}*U.cm,
        inclH2=true,
        inclH2O=false,
        inclReact=false,
        inclYFaces=false,
        inclZFaces=false,
        inclXFaces=true,
        inclVelX=true,
        H2(
          xNegative(
            slipY=false,
            slipZ=false,
            thermoOpt=ThermoOpt.OpenDiabatic),
          xPositive(
            slipY=false,
            slipZ=false,
            final thermoOpt=ThermoOpt.ClosedAdiabatic),
          setVelX=true))
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));
      inner Defaults defaults(analysis=true, T=293.15*U.K)
        annotation (Placement(transformation(extent={{-60,70},{-40,90}})));

    protected
      Connectors.FaceBusInternal xNegative "Positive face along the x axis"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    equation

      connect(xNegative.H2.material, volume.face.material) annotation (Line(
          points={{5.55112e-16,5.55112e-16},{-6,5.55112e-16},{-6,0},{-10,0},{-10,
              6.10623e-16},{-20,6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(xNegative.H2.thermal, volume.face.thermal) annotation (Line(
          points={{5.55112e-16,5.55112e-16},{-6,5.55112e-16},{-6,0},{-10,0},{-10,
              6.10623e-16},{-20,6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(xNegative, subregion.xNegative) annotation (Line(
          points={{5.55112e-16,5.55112e-16},{6,5.55112e-16},{6,0},{10,0},{10,
              6.10623e-16},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.FluidAdapt.mos"));
    end AdaptFluid2;
  end Examples;

  package Adapters "Adapters to Package Modelica"
    extends Modelica.Icons.Package;

    model AdaptBusH2
      "<html>Adapter for H<sub>2</sub> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
      extends BaseClasses.PartialAdaptBus(redeclare replaceable package Medium
          = Modelica.Media.IdealGases.SingleGases.H2 (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false));

    equation
      connect(fluidAdapt.face.material, face.H2.material) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));

      connect(fluidAdapt.face.thermal, face.H2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));

      annotation (Diagram(graphics));
    end AdaptBusH2;

    model AdaptBusH2O
      "<html>Adapter for H<sub>2</sub>O between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
      extends BaseClasses.PartialAdaptBus(redeclare replaceable package Medium
          = Modelica.Media.IdealGases.SingleGases.H2O);

    equation
      connect(fluidAdapt.face.material, face.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));

      connect(fluidAdapt.face.thermal, face.H2O.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));
    end AdaptBusH2O;

    model AdaptBusN2
      "<html>Adapter for N<sub>2</sub> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
      extends BaseClasses.PartialAdaptBus(redeclare replaceable package Medium
          = Modelica.Media.IdealGases.SingleGases.N2);
    equation
      connect(fluidAdapt.face.material, face.N2.material) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));

      connect(fluidAdapt.face.thermal, face.N2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));
    end AdaptBusN2;

    model AdaptBusO2
      "<html>Adapter for O<sub>2</sub> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
      extends BaseClasses.PartialAdaptBus(redeclare replaceable package Medium
          = Modelica.Media.IdealGases.SingleGases.O2);

    equation
      connect(fluidAdapt.face.material, face.O2.material) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));

      connect(fluidAdapt.face.thermal, face.O2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-10,3},{-10,3}}));
    end AdaptBusO2;

    model AdaptFluid
      "<html>Fluid adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
      extends FCSys.BaseClasses.Icons.Names.Top3;

      // TODO:  Fix this; there may be issues with the energy/entropy rate balance.

      parameter Q.Area A=1*U.cm^2 "Area of the connection surface";
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium (
            final nXi=0) "Medium model" annotation (choicesAllMatching=true);

      FCSys.Connectors.FaceGeneric face(
        final thermoOpt=ThermoOpt.OpenDiabatic,
        final slip1=false,
        final slip2=false)
        "Connector for material and entropy of a single species" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
          Medium = Medium) "Modelica fluid port" annotation (Placement(
            transformation(extent={{90,-10},{110,10}}), iconTransformation(
              extent={{90,-10},{110,10}})));

      Medium.BaseProperties medium "Base properties of the fluid";
      Q.NumberAbsolute s "Specific entropy";

    equation
      // Thermodynamic state and properties
      medium.p = fluidPort.p;
      medium.T = face.thermal.T/U.K;
      medium.Xi = ones(Medium.nXi)/Medium.nXi;
      s = Medium.specificEntropy(medium.state)*medium.MM*U.J/(U.mol*U.K);

      // Efforts and streams
      //  face.material.mu = actualStream(fluidPort.h_outflow)*medium.MM*U.J/U.mol -  face.thermal.T*s;
      //  face.material.mu = inStream(fluidPort.h_outflow)*medium.MM*U.J/U.mol - face.thermal.T*s;
      face.material.mu = medium.h*medium.MM*U.J/U.mol - face.thermal.T*s;
      fluidPort.h_outflow = medium.h;

      // Rate balances (no storage)
      0 = (medium.MM*U.kg/U.mol)*face.material.Ndot + fluidPort.m_flow*U.kg/U.s
        "Mass";
      face.thermal.Qdot = s*face.material.Ndot
        "No thermal conduction--advection only";
      // The rest of the energy balance cancels.

      annotation (
        Documentation(info="<html><p>Note that transverse momentum is not included.</p>
  </html>"),
        Diagram(graphics),
        Icon(graphics={Line(
                  points={{0,40},{0,-40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-100,0}},
                  color={127,127,127},
                  smooth=Smooth.None),Line(
                  points={{0,0},{100,0}},
                  color={0,127,255},
                  smooth=Smooth.None)}));
    end AdaptFluid;

    model 'AdaptSubregione-'
      "<html>Electrical adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> with bus connector</html>"
      import FCSys.BaseClasses.Axis;
      extends FCSys.BaseClasses.Icons.Names.Top3;

      parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

      FCSys.Connectors.FaceBus face
        "Connector for material, linear momentum, and entropy" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {90,30},{110,50}}),iconTransformation(extent={{90,30},{110,50}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port
        "Modelica heat port" annotation (Placement(transformation(extent={{90,-50},
                {110,-30}}), iconTransformation(extent={{90,-50},{110,-30}})));
      FCSys.BCs.Adapters.'AdaptBuse-' electAdapt
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      connect(electAdapt.face, face.graphite) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(electAdapt.pin, pin) annotation (Line(
          points={{10,4},{80,4},{80,40},{100,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electAdapt.port, port) annotation (Line(
          points={{10,-4},{80,-4},{80,-40},{100,-40}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (Diagram(graphics), Icon(graphics={Line(
                  points={{0,40},{0,-40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-100,0}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,40},{100,40}},
                  color={0,0,255},
                  smooth=Smooth.None),Line(
                  points={{0,-40},{100,-40}},
                  color={191,0,0},
                  smooth=Smooth.None)}));
    end 'AdaptSubregione-';

    model 'AdaptBuse-'
      "<html>Electrical adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> with bus connector</html>"
      import FCSys.BaseClasses.Axis;
      extends FCSys.BaseClasses.Icons.Names.Top3;

      parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

      FCSys.Connectors.FaceBus face
        "Connector for material, linear momentum, and entropy" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {90,30},{110,50}}),iconTransformation(extent={{90,30},{110,50}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port
        "Modelica heat port" annotation (Placement(transformation(extent={{90,-50},
                {110,-30}}), iconTransformation(extent={{90,-50},{110,-30}})));
      FCSys.BCs.Adapters.'Adapte-' electAdapt
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      connect(electAdapt.face.material, face.'e-'.material) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(electAdapt.face.thermal, face.'e-'.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(electAdapt.pin, pin) annotation (Line(
          points={{10,4},{80,4},{80,40},{100,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(electAdapt.port, port) annotation (Line(
          points={{10,-4},{80,-4},{80,-40},{100,-40}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (Diagram(graphics), Icon(graphics={Line(
                  points={{0,40},{0,-40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-100,0}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,40},{100,40}},
                  color={0,0,255},
                  smooth=Smooth.None),Line(
                  points={{0,-40},{100,-40}},
                  color={191,0,0},
                  smooth=Smooth.None)}));
    end 'AdaptBuse-';

    model 'Adapte-'
      "<html>Electrical adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (<a href=\"modelica://Modelica.Electrical.Analog\">Electrical.Analog</a> and <a href=\"modelica://Modelica.Thermal.HeatTransfer\">Thermal.HeatTransfer</a>)</html>"
      import Data = FCSys.Characteristics.'e-'.Gas;

      extends FCSys.BaseClasses.Icons.Names.Top3;

      Connectors.FaceGeneric face(
        final thermoOpt=ThermoOpt.OpenDiabatic,
        final slip1=false,
        final slip2=false)
        "Connector for material, linear momentum, and entropy of a single species"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {90,30},{110,50}}),iconTransformation(extent={{90,30},{110,50}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port
        "Modelica heat port" annotation (Placement(transformation(extent={{90,-50},
                {110,-30}}), iconTransformation(extent={{90,-50},{110,-30}})));

    equation
      // Equal efforts
      face.material.mu = pin.v*U.V "Electrochemical potential";
      face.thermal.T = port.T*U.K "Temperature";

      // Conservation (no storage)
      0 = face.material.Ndot - pin.i*U.A "Material";
      // In FCSys, current is material (e.g., electron) flow, not charge flow;
      // therefore, the flow rate is negated.
      0 = face.thermal.Qdot + port.Q_flow*U.W "Energy";
      // There is no electrical work since electrons are not stored and there
      // is no potential difference.

      annotation (Diagram(graphics), Icon(graphics={Line(
                  points={{0,40},{0,-40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-100,0}},
                  color={127,127,127},
                  smooth=Smooth.None),Line(
                  points={{0,40},{100,40}},
                  color={0,0,255},
                  smooth=Smooth.None),Line(
                  points={{0,-40},{100,-40}},
                  color={191,0,0},
                  smooth=Smooth.None)}));
    end 'Adapte-';

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      model PartialAdaptBus
        "<html>Partial model for a fluid adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> with bus connector</html>"
        extends FCSys.BaseClasses.Icons.Names.Top3;

        parameter Q.Area A=1*U.cm^2 "Area of the connection surface";
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
          constrainedby Modelica.Media.Interfaces.PartialMedium "Medium model"
          annotation (choicesAllMatching=true);

        FCSys.Connectors.FaceBus face
          "Connector for material, linear momentum, and entropy" annotation (
            Placement(transformation(extent={{-110,-10},{-90,10}}),
              iconTransformation(extent={{-110,-10},{-90,10}})));
        FCSys.BCs.Adapters.AdaptFluid fluidAdapt(final A=A, redeclare final
            package Medium = Medium)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{90,-10},{110,10}})));

      equation
        connect(fluidAdapt.fluidPort, fluidPort) annotation (Line(
            points={{10,6.10623e-16},{80,6.10623e-16},{80,5.55112e-16},{100,
                5.55112e-16}},
            color={0,127,255},
            smooth=Smooth.None));

        annotation (Diagram(graphics), Icon(graphics={Line(
                      points={{0,40},{0,-40}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{-100,0}},
                      color={127,127,127},
                      smooth=Smooth.None),Line(
                      points={{0,0},{100,0}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end PartialAdaptBus;
    end BaseClasses;
  end Adapters;

  package TestStands "Test stands"
    extends Modelica.Icons.Package;
    model TestProfile "Test profile"
      extends Modelica.Icons.Example;
      extends FCSys.BCs.TestStands.BaseClasses.PartialTestStandNoIO;
      /*(
    anEndBC(each graphite(inclC=true, 'incle-'=true)),
    anSourceBC(each gas(inclH2=true, inclH2O=true)),
    anSinkBC(each gas(inclH2=true, inclH2O=true)),
    caEndBC(each graphite('incle-'=true)),
    caSourceBC(each gas(
        inclH2O=true,
        inclN2=true,
        inclO2=true)),
    caSinkBC(each graphite(inclC=true, 'incle-'=true)));
*/

      annotation (
        defaultComponentName="testStand",
        defaultComponentPrefixes="replaceable",
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}),graphics));
    end TestProfile;

    model ReplayData
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
        fileName="FCSys/tests/LOOCV/data.mat",
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
            origin={-148,0})));
      Modelica.Blocks.Math.Gain unit2(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-130,-30})));
      Modelica.Blocks.Math.Gain unit3(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-110,0})));
      FCSys.BCs.BaseClasses.RealFunction unit4(y=unit4.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit5(y=unit5.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,0})));
      FCSys.BCs.BaseClasses.RealFunction unit6(y=unit6.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit7(y=unit7.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,0})));
      Modelica.Blocks.Math.Gain unit8(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,-30})));
      Modelica.Blocks.Math.Gain unit9(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,0})));
      FCSys.BCs.BaseClasses.RealFunction unit10(y=(unit10.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit11(y=(unit11.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,0})));
      FCSys.BCs.BaseClasses.RealFunction unit12(y=(unit12.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={70,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit13(y=(unit13.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,0})));
      FCSys.BCs.BaseClasses.RealFunction unit14(y=(unit14.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={110,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit15(y=(unit15.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={130,0})));
      Modelica.Blocks.Math.Gain unit16(k=U.A) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={150,-30})));

      Connectors.RealOutputInternal Deltamu(final unit="l2.m/(N.T2)")
        "CVM Cell 1 Voltage" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-148,-60})));
      Connectors.RealOutputInternal RHAnFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Anode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-130,-60})));
      Connectors.RealOutputInternal RHCaFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Cathode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-110,-60})));
      Connectors.RealOutputInternal p_anFPNegY(final unit="m/(l.T2)")
        "Pressure anode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,-60})));
      Connectors.RealOutputInternal p_anFPPosY(final unit="m/(l.T2)",final min=
            0) "Pressure anode outlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,-60})));
      Connectors.RealOutputInternal p_caFPNegY(final unit="m/(l.T2)")
        "Pressure cathode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,-60})));
      Connectors.RealOutputInternal p_caFPPosY(final unit="m/(l.T2)",final min=
            0) "Pressure anode outlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,-60})));
      Connectors.RealOutputInternal Vdot_anFPNegY_H2(final unit="l3/T")
        "Flow anode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,-60})));
      Connectors.RealOutputInternal Vdot_caFPNegY_air(final unit="l3/T")
        "Flow cathode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,-60})));
      Connectors.RealOutputInternal T_anFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-60})));
      Connectors.RealOutputInternal T_anFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,-60})));
      Connectors.RealOutputInternal T_caFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={70,-60})));
      Connectors.RealOutputInternal T_caFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,-60})));
      Connectors.RealOutputInternal T_anFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate anode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={110,-60})));
      Connectors.RealOutputInternal T_caFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate cathode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={130,-60})));
      Connectors.RealOutputInternal J(final unit="N/T") "Measured load"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={150,-60})));

    equation
      //  Terminate as desired
      if terminateMaxTime then
        when time > combiTimeTable.t_max then
          terminate("The end of the data has been reached.");
        end when;
      end if;

      // Connections from source to unit conversion
      connect(unit1.u, combiTimeTable.y[1]) annotation (Line(
          points={{-148,12},{-148,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit2.u, combiTimeTable.y[2]) annotation (Line(
          points={{-130,-18},{-130,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit3.u, combiTimeTable.y[3]) annotation (Line(
          points={{-110,12},{-110,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit4.u, combiTimeTable.y[4]) annotation (Line(
          points={{-90,-20},{-90,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit5.u, combiTimeTable.y[5]) annotation (Line(
          points={{-70,10},{-70,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit6.u, combiTimeTable.y[6]) annotation (Line(
          points={{-50,-20},{-50,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit7.u, combiTimeTable.y[7]) annotation (Line(
          points={{-30,10},{-30,50},{-1.44329e-15,50},{-1.44329e-15,69}},
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
          points={{30,-20},{30,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit11.u, combiTimeTable.y[13]) annotation (Line(
          points={{50,10},{50,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit12.u, combiTimeTable.y[14]) annotation (Line(
          points={{70,-20},{70,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit13.u, combiTimeTable.y[15]) annotation (Line(
          points={{90,10},{90,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit14.u, combiTimeTable.y[16]) annotation (Line(
          points={{110,-20},{110,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit15.u, combiTimeTable.y[17]) annotation (Line(
          points={{130,10},{130,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit16.u, combiTimeTable.y[18]) annotation (Line(
          points={{150,-18},{150,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from unit conversion to internal outputs
      connect(Deltamu, unit1.y) annotation (Line(
          points={{-148,-60},{-148,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHAnFPNegX, unit2.y) annotation (Line(
          points={{-130,-60},{-130,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHCaFPNegX, unit3.y) annotation (Line(
          points={{-110,-60},{-110,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPNegY, unit4.y) annotation (Line(
          points={{-90,-60},{-90,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPPosY, unit5.y) annotation (Line(
          points={{-70,-60},{-70,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPNegY, unit6.y) annotation (Line(
          points={{-50,-60},{-50,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPPosY, unit7.y) annotation (Line(
          points={{-30,-60},{-30,-35},{-30,-10},{-30,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_anFPNegY_H2, unit8.y) annotation (Line(
          points={{-10,-60},{-10,-50.5},{-10,-41},{-10,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_caFPNegY_air, unit9.y) annotation (Line(
          points={{10,-60},{10,-40},{14,-40},{14,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, unit10.y) annotation (Line(
          points={{30,-60},{30,50},{30,-40},{30,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, unit11.y) annotation (Line(
          points={{50,-60},{50,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, unit12.y) annotation (Line(
          points={{70,-60},{70,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, unit13.y) annotation (Line(
          points={{90,-60},{90,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, unit14.y) annotation (Line(
          points={{110,-60},{110,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPX, unit15.y) annotation (Line(
          points={{130,-60},{130,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(J, unit16.y) annotation (Line(
          points={{150,-60},{150,-41}},
          color={0,0,127},
          smooth=Smooth.None));

      // Summations
      connect(sumAnMFC.y, unit8.u) annotation (Line(
          points={{-14,19},{-14,-10},{-10,-10},{-10,-18}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sumCaMFC.y, unit9.u) annotation (Line(
          points={{14,19},{14,17.25},{14,17.25},{14,15.5},{14,12},{14,12}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from internal outputs to public output
      connect(Deltamu, y.Deltamu) annotation (Line(
          points={{-148,-60},{-148,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHAnFPNegX, y.RHAnFPNegX) annotation (Line(
          points={{-130,-60},{-130,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHCaFPNegX, y.RHCaFPNegX) annotation (Line(
          points={{-110,-60},{-110,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPNegY, y.p_anFPNegY) annotation (Line(
          points={{-90,-60},{-90,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPPosY, y.p_anFPPosY) annotation (Line(
          points={{-70,-60},{-70,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPNegY, y.p_caFPNegY) annotation (Line(
          points={{-50,-60},{-50,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPPosY, y.p_caFPPosY) annotation (Line(
          points={{-30,-60},{-30,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_anFPNegY_H2, y.Vdot_anFPNegY_H2) annotation (Line(
          points={{-10,-60},{-10,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_caFPNegY_air, y.Vdot_caFPNegY_air) annotation (Line(
          points={{10,-60},{10,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, y.T_anFPNegY) annotation (Line(
          points={{30,-60},{30,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, y.T_anFPPosY) annotation (Line(
          points={{50,-60},{50,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, y.T_caFPNegY) annotation (Line(
          points={{70,-60},{70,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, y.T_caFPPosY) annotation (Line(
          points={{90,-60},{90,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, y.T_anFPX) annotation (Line(
          points={{110,-60},{110,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T_caFPX, y.T_caFPX) annotation (Line(
          points={{130,-60},{130,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(J, y.J) annotation (Line(
          points={{150,-60},{150,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (
        Commands(file="tests/LOOCV/LOOCV.mos"
            "Perform leave-one-out cross validation on the cell model"),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-180,-100},
                {180,100}}),graphics),
        Icon(graphics),
        experiment(StopTime=15481, Algorithm="Euler"),
        experimentSetupOutput);
    end ReplayData;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialTestStand "Partial test stand"
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
              origin={-160,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-160,0})));
        FCSys.Connectors.FaceBus caEnd[n_y, n_z] "Cathode end plate"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={160,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={162,0})));
        FCSys.Connectors.FaceBus anSource[n_x_an, n_z] "Anode source"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-160}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,-160})));
        FCSys.Connectors.FaceBus anSink[n_x_an, n_z] "Anode sink" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,160}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,160})));
        FCSys.Connectors.FaceBus caSource[n_x_ca, n_z] "Cathode source"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,160})));
        FCSys.Connectors.FaceBus caSink[n_x_ca, n_z] "Cathode sink" annotation
          (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,160}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,-160})));

        FCSys.BCs.Face.Subregion0Current anEndBC[n_y, n_z](each final axis=
              FCSys.BaseClasses.Axis.x, each graphite(inclC=true, 'incle-'=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-136,0})));
        FCSys.BCs.Face.Subregion0Current caEndBC[n_y, n_z](each final axis=
              FCSys.BaseClasses.Axis.x, each graphite(inclC=true, 'incle-'=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={136,0})));
        FCSys.BCs.Face.Subregion0Current anSourceBC[n_x_an, n_z](each final
            axis=FCSys.BaseClasses.Axis.y, each gas(inclH2=true, inclH2O=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-136})));
        FCSys.BCs.Face.Subregion0Current anSinkBC[n_x_an, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(inclH2=true, inclH2O=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-40,136})));
        FCSys.BCs.Face.Subregion0Current caSourceBC[n_x_ca, n_z](each final
            axis=FCSys.BaseClasses.Axis.y, each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-136})));
        FCSys.BCs.Face.Subregion0Current caSinkBC[n_x_ca, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(
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
              origin={-160,160})));
        Connectors.RealOutputBus y[n_y, n_z] if inclIO annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={160,-160})));
        replaceable BCs.FaceDifferential.Subregion current[n_y, n_z](each
            final axis=FCSys.BaseClasses.Axis.x, graphite(
            inclC=true,
            C(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic),

            'incle-'=true,
            'e-'(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic)))
          if inclIO constrainedby BCs.FaceDifferential.Subregion(graphite(
            inclC=true,
            C(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.ClosedAdiabatic),

            'incle-'=true,
            'e-'(thermoOpt=FCSys.Connectors.BaseClasses.ThermoOpt.OpenDiabatic)))
          annotation (Placement(transformation(extent={{-140,20},{-120,40}})));

        replaceable Sensors.FaceDifferential.Subregion voltage[n_y, n_z](each
            final axis=FCSys.BaseClasses.Axis.x) if inclIO
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
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},
                  {160,160}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{
                  160,160}}), graphics={Rectangle(
                      extent={{-160,160},{160,-160}},
                      lineColor={191,191,191},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Backward),Rectangle(extent={{-160,
                160},{160,-160}}, lineColor={0,0,0})}));
      end PartialTestStand;

      partial model PartialTestStandNoIO "Partial test stand"
        extends FCSys.BaseClasses.Icons.Names.Top9;
        final parameter Integer n_x_an=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x an</sub>)</html>";
        final parameter Integer n_x_ca=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x ca</sub>)</html>";
        final parameter Integer n_y=1
          "<html>Number of subregions along the channel (<i>n</i><sub>y</sub>)</html>";
        final parameter Integer n_z=1
          "<html>Number of subregions across the channel (<i>n</i><sub>z</sub>)</html>";

        FCSys.BCs.Face.Subregion0Current anEnd[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite(
            inclC=true,
            'incle-'=true,
            'e-'(redeclare FCSys.BCs.Face.Species.Material.Current materialBC,
                redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=1*U.A,
                  duration=50)))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-30,0})));
        FCSys.BCs.Face.Subregion0Current caEnd[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite(
            inclC=true,
            'incle-'=true,
            'e-'(redeclare FCSys.BCs.Face.Species.Material.Current materialBC,
                redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=-1*U.A,
                  duration=50)))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,0})));
        FCSys.BCs.Face.Subregion0Current anSource[n_x_an, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(inclH2=true, inclH2O=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-30})));
        FCSys.BCs.Face.Subregion0Current anSink[n_x_an, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(inclH2=true, inclH2O=true))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,30})));
        FCSys.BCs.Face.Subregion0Current caSource[n_x_ca, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,-30})));
        FCSys.BCs.Face.Subregion0Current caSink[n_x_ca, n_z](each final axis=
              FCSys.BaseClasses.Axis.y, each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={20,30})));

        inner Defaults defaults
          annotation (Placement(transformation(extent={{50,20},{70,40}})));
      equation

        annotation (
          defaultComponentName="testStand",
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{
                  160,160}}), graphics));
      end PartialTestStandNoIO;
    end BaseClasses;
  end TestStands;

  package Chemical
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> and <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connectors</html>"
    extends Modelica.Icons.Package;

    package Phases
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;
      model Phase "BC for a phase with all species conditionally included"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species C(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.C.Graphite) if inclC "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));

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
        FCSys.BCs.Chemical.Species.Species C19HF37O5S(final n_vel=n_vel,
            redeclare package Data = FCSys.Characteristics.C19HF37O5S.Solid)
          if inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species 'e-'(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.'e-'.Graphite) if 'incle-'
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species H2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2.Gas) if inclH2 "Model"
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
        FCSys.BCs.Chemical.Species.Species H2O(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2O.Gas) if inclH2O "Model"
          annotation (Dialog(
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
        FCSys.BCs.Chemical.Species.Species 'H+'(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.'H+'.Solid) if 'inclH+'
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species N2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.N2.Gas) if inclN2 "Model"
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
        FCSys.BCs.Chemical.Species.Species O2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.O2.Gas) if inclO2 "Model"
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

        // C
        connect(C.chemical, chemical.C) annotation (Line(
            points={{-5.08852e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
                4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
            points={{-5.08852e-16,-4},{5.55112e-16,-40}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
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

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalBC");
      end Phase;

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
        FCSys.BCs.Chemical.Species.Species H2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2.Gas) if inclH2 "Model"
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
        FCSys.BCs.Chemical.Species.Species H2O(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2O.Gas) if inclH2O "Model"
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
        FCSys.BCs.Chemical.Species.Species N2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.N2.Gas) if inclN2 "Model"
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
        FCSys.BCs.Chemical.Species.Species O2(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.O2.Gas) if inclO2 "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
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

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalBC");
      end Gas;

      model Graphite "BC for graphite"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species C(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.C.Graphite) if inclC "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Chemical.Species.Species 'e-'(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.'e-'.Graphite) if 'incle-'
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // Note:  It would be helpful if Modelica allowed elements of expandable
        // connectors to be named by the contents of a string variable and the
        // name of an instance of a model was accessible through a string (like
        // %name is expanded to be the name of the instance of the model).  Then,
        // the connection equations that follow could be generic.

        // C
        connect(C.chemical, chemical.C) annotation (Line(
            points={{-5.08852e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
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

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalBC");
      end Graphite;

      model Ionomer "BC for ionomer"

        extends BaseClasses.NullPhase;

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
        FCSys.BCs.Chemical.Species.Species C19HF37O5S(final n_vel=n_vel,
            redeclare package Data = FCSys.Characteristics.C19HF37O5S.Solid)
          if inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
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
        FCSys.BCs.Chemical.Species.Species H2O(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2O.Gas) if inclH2O "Model"
          annotation (Dialog(
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
        FCSys.BCs.Chemical.Species.Species 'H+'(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.'H+'.Solid) if 'inclH+'
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
            points={{-5.08852e-16,-4},{5.55112e-16,-40}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
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

        annotation (
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2),
          Placement(transformation(extent={{-10,-10},{10,10}})),
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalBC",
          Diagram(graphics));
      end Ionomer;

      model Liquid "BC for liquid"

        extends BaseClasses.NullPhase;

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
        FCSys.BCs.Chemical.Species.Species C19HF37O5S(final n_vel=n_vel,
            redeclare package Data = FCSys.Characteristics.C19HF37O5S.Solid)
          if inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
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
        FCSys.BCs.Chemical.Species.Species H2O(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.H2O.Gas) if inclH2O "Model"
          annotation (Dialog(
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
        FCSys.BCs.Chemical.Species.Species 'H+'(final n_vel=n_vel, redeclare
            package Data = FCSys.Characteristics.'H+'.Solid) if 'inclH+'
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
            points={{-5.08852e-16,-4},{5.55112e-16,-40}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
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

        annotation (
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2),
          Placement(transformation(extent={{-10,-10},{10,10}})),
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalBC",
          Diagram(graphics));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        model NullPhase "Empty BC for a phase (no species)"
          extends FCSys.BaseClasses.Icons.BCs.Single;
          parameter Integer n_vel(
            final min=1,
            final max=3) = 1
            "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
            annotation (HideResult=true);

          FCSys.Connectors.ChemicalBus chemical
            "Multi-species connector for material"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,40})));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalBC",
            Diagram(graphics));
        end NullPhase;
      end BaseClasses;
    end Phases;

    package Species
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends Modelica.Icons.Package;

      model Species
        "<html>BCs for the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

        import FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial;
        import FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMomentum;
        import FCSys.BCs.Chemical.Species.BaseClasses.BCTypeEnthalpy;
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Integer n_vel(
          final min=1,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);
        replaceable package Data =
            FCSys.Characteristics.BaseClasses.Characteristic constrainedby
          FCSys.Characteristics.BaseClasses.Characteristic
          "Characteristic data of the species" annotation (
          Dialog(group="Material properties"),
          __Dymola_choicesAllMatching=true,
          Placement(transformation(extent={{-60,40},{-40,60}}),
              iconTransformation(extent={{-10,90},{10,110}})));

        parameter BCTypeMaterial materialBC=BCTypeMaterial.PotentialElectrochemicalPerTemperature
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
        replaceable Modelica.Blocks.Sources.Constant materialSpec(k(start=
                Data.g(300*U.K)/(300*U.K))) if internalMaterial constrainedby
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

        // 1st component of velocity
        parameter BCTypeMomentum vel1BC=BCTypeMomentum.Velocity "Type of BC"
          annotation (Dialog(
            group="1st component of velocity",
            enable=n_vel > 0,
            __Dymola_descriptionLabel=true));
        parameter Boolean internalVel1=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="1st component of velocity",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant vel1Spec(k(start=0)) if
          internalVel1 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st component of velocity",
            __Dymola_descriptionLabel=true,
            enable=internalVel1),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-50,10})));

        // 2nd component of velocity
        parameter BCTypeMomentum vel2BC=BCTypeMomentum.Velocity "Type of BC"
          annotation (Dialog(
            group="2nd component of velocity",
            enable=n_vel > 1,
            __Dymola_descriptionLabel=true));
        parameter Boolean internalVel2=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="2nd component of velocity",
            enable=n_vel > 1,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant vel2Spec(k(start=0)) if
          internalVel2 and n_vel > 1 constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd component of velocity",
            enable=n_vel > 1,
            __Dymola_descriptionLabel=true,
            enable=internalVel2),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,10})));

        // 3rd component of velocity
        parameter BCTypeMomentum vel3BC=BCTypeMomentum.Velocity "Type of BC"
          annotation (Dialog(
            group="3rd component of velocity",
            enable=n_vel > 2,
            __Dymola_descriptionLabel=true));
        parameter Boolean internalVel3=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="3rd component of velocity",
            enable=n_vel > 2,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant vel3Spec(k(start=0)) if
          internalVel3 and n_vel > 2 constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="3rd component of velocity",
            enable=n_vel > 2,
            __Dymola_descriptionLabel=true,
            enable=internalVel3),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,10})));

        // Enthalpy
        parameter BCTypeMomentum enthalpyBC=BCTypeEnthalpy.EnthalpyMassic
          "Type of BC"
          annotation (Dialog(group="Enthalpy", __Dymola_descriptionLabel=true));
        parameter Boolean internalEnth=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Enthalpy",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant enthSpec(k(start=Data.h()))
          if internalEnth constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Enthalpy",
            __Dymola_descriptionLabel=true,
            enable=internalEnth),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,10})));

        Connectors.ChemicalOutput chemical(
          final n_vel=n_vel,
          final formula=Data.formula,
          final m=Data.m)
          "Single-species connector for material, with advection of linear momentum and enthalpy"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
        Connectors.RealInputBus u "Bus of inputs to specify conditions"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        Connectors.RealInputInternal u_N(final unit=if materialBC ==
              BCTypeMaterial.PotentialElectrochemicalPerTemperature then "1"
               else "N/T") if not internalMaterial "Material signal"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-70,30}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_Phi_1(final unit="l/T") if not
          internalVel1 and n_vel > 0 "Signal for the 1st component of velocity"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,30}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_Phi_2(final unit="l/T") if not
          internalVel2 and n_vel > 1 "Signal for the 2nd component of velocity"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_Phi3(final unit="l/T") if not
          internalVel3 and n_vel > 2 "Signal for the 3rd component of velocity"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={50,30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_H(final unit="m.l2/(N.T2)") if not
          internalEnth "Signal for enthalpy" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={90,30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_N_int "Internal material signal"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-70,-20}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_Phi_int[n_vel]
          "Internal signal for linear momentum" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,-20}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        Connectors.RealInputInternal u_H_int "Internal signal for enthalpy"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={90,-20}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      equation
        // Material
        if materialBC == BCTypeMaterial.PotentialElectrochemicalPerTemperature
             then
          chemical.muPerT = u_N_int;
        else
          chemical.Ndot = u_N_int;
        end if;

        // 1st component of velocity
        if n_vel > 0 then
          //  if bCTypeMomentum1 == BCTypeMomentum.Velocity then
          chemical.phi[1] = u_Phi_int[1];
          //  end if;
        end if;

        // 2nd component of velocity
        if n_vel > 1 then
          //  if bCTypeMomentum2 == BCTypeMomentum.Velocity then
          chemical.phi[2] = u_Phi_int[2];
          //  end if;
        end if;

        // 3rd component of velocity
        if n_vel > 2 then
          //  if bCTypeMomentum3 == BCTypeMomentum.Velocity then
          chemical.phi[3] = u_Phi_int[3];
          //  end if;
        end if;

        // Enthalpy
        //  if enthalpyBC== BCTypeEnthalpy.EnthalpyMassic then
        chemical.hbar = u_H_int;
        //  end if;

        connect(u.N, u_N) annotation (Line(
            points={{5.55112e-16,60},{5.55112e-16,40},{-70,40},{-70,30}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(u.Phi_1, u_Phi_1) annotation (Line(
            points={{5.55112e-16,60},{5.55112e-16,40},{-30,40},{-30,30}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(u.Phi_2, u_Phi_2) annotation (Line(
            points={{5.55112e-16,60},{5.55112e-16,40},{10,40},{10,30}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(u.Phi3, u_Phi3) annotation (Line(
            points={{5.55112e-16,60},{5.55112e-16,40},{50,40},{50,30}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(u.H, u_H) annotation (Line(
            points={{5.55112e-16,60},{5.55112e-16,40},{90,40},{90,30}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(u_N_int, u_N) annotation (Line(
            points={{-70,-20},{-70,30}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_1, u_Phi_int[1]) annotation (Line(
            points={{-30,30},{-30,-8},{10,-8},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_2, u_Phi_int[2]) annotation (Line(
            points={{10,30},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi3, u_Phi_int[3]) annotation (Line(
            points={{50,30},{50,-8},{10,-8},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_H, u_H_int) annotation (Line(
            points={{90,30},{90,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(materialSpec.y, u_N_int) annotation (Line(
            points={{-90,-1},{-90,-8},{-70,-8},{-70,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(vel1Spec.y, u_Phi_int[1]) annotation (Line(
            points={{-50,-1},{-50,-8},{10,-8},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(vel2Spec.y, u_Phi_int[2]) annotation (Line(
            points={{30,-1},{30,-8},{10,-8},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(vel3Spec.y, u_Phi_int[3]) annotation (Line(
            points={{-10,-1},{-10,-8},{10,-8},{10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(enthSpec.y, u_H_int) annotation (Line(
            points={{70,-1},{70,-8},{90,-8},{90,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (
          defaultComponentName="speciesChemicalBC",
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,-100},
                  {120,100}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-120,-100},{
                  120,100}}), graphics));
      end Species;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        type BCTypeMaterial = enumeration(
            PotentialElectrochemicalPerTemperature
              "Prescribed quotient of electrochemical potential and temperature",

            Current "Prescribed current") "Types of BCs for material";

        type BCTypeMomentum = enumeration(
            Velocity "Prescribed velocity") "Types of BCs for linear momentum";
        type BCTypeEnthalpy = enumeration(
            EnthalpyMassic "Prescribed massic enthalpy")
          "Types of BCs for enthalpy";
      end BaseClasses;
    end Species;
    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the boundary conditions must be as well.  A
<a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a>
connector
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.BCs.Chemical.Species.Species\">Species
boundary condition</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.BCs.Chemical.Phases\">Phase
boundary condition</a> models.
</p></html>"));
  end Chemical;

  package InertAmagat
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector</html>"
    extends Modelica.Icons.Package;

    model Phase
      "<html>BC for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> model</html>"
      extends FCSys.BaseClasses.Icons.BCs.Single;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>";

      // Additivity of pressure
      replaceable FCSys.BCs.InertAmagat.Volume.Volume volBC(final n_vel=n_vel)
        constrainedby FCSys.BCs.InertAmagat.Volume.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Volume", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-90,-26},{-70,-6}})));
      parameter Boolean internalVol=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Volume",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant volSpec(k(start=1*U.cm^3))
        if internalVol constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Volume",
          __Dymola_descriptionLabel=true,
          enable=internalVol),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,20})));

      // 1st component of velocity
      replaceable FCSys.BCs.InertAmagat.Momentum.Force lin1BC(final n_vel=n_vel)
        constrainedby FCSys.BCs.InertAmagat.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="1st axis",
          enable=n_vel > 0,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-50,-26},{-30,-6}})));
      parameter Boolean internalLin1=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="1st axis",
          enable=n_vel > 0,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin1Spec(k(start=0)) if
        internalLin1 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal signal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="1st axis",
          enable=n_vel > 0,
          __Dymola_descriptionLabel=true,
          enable=internalLin1),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,20})));

      // 2nd component of velocity
      replaceable FCSys.BCs.InertAmagat.Momentum.Force lin2BC(final n_vel=n_vel)
        if n_vel > 1 constrainedby
        FCSys.BCs.InertAmagat.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-26},{10,-6}})));
      parameter Boolean internalLin2=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin2Spec(k(start=0)) if
        internalLin2 and n_vel > 1 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1 and internalLin2,
          __Dymola_descriptionLabel=true,
          enable=internalLin2),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,20})));

      // 3rd component of velocity
      replaceable FCSys.BCs.InertAmagat.Momentum.Force lin3BC(final n_vel=n_vel)
        if n_vel > 2 constrainedby
        FCSys.BCs.InertAmagat.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{30,-26},{50,-6}})));
      parameter Boolean internalLin3=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin3Spec(k(start=0)) if
        internalLin3 and n_vel > 2 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2 and internalLin3,
          __Dymola_descriptionLabel=true,
          enable=internalLin3),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,20})));

      // Heat
      replaceable FCSys.BCs.InertAmagat.Heat.HeatFlowRate heatBC(final n_vel=
            n_vel) constrainedby
        FCSys.BCs.InertAmagat.Heat.BaseClasses.PartialBC "Type of condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Heat", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{70,-26},{90,-6}})));
      parameter Boolean internalHeat=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Heat",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant heatSpec(k(start=0)) if
        internalHeat constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Heat",
          __Dymola_descriptionLabel=true,
          enable=internalHeat),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={70,20})));

      Connectors.RealInputBus u "Input bus for external signal sources"
        annotation (HideResult=not (internalVol or internalLin1 or internalLin2
             or internalLin3 or internalHeat), Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,70}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

      Connectors.InertAmagat inert(final n_vel=n_vel)
        "Single-species connector for linear momentum and heat, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
    protected
      Connectors.RealInputInternal u_N if not internalVol
        "External material signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,40})));
      Connectors.RealInputInternal u_Phi_1 if not internalLin1
        "External signal for 1st axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Phi_2 if not internalLin2 and n_vel > 1
        "External signal for 2nd axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Phi3 if not internalLin3 and n_vel > 2
        "External signal for 3rd axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Q if not internalHeat
        "External signal for heat" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,40}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,40})));

    equation
      // Material
      connect(volBC.inert, inert) annotation (Line(
          points={{-80,-20},{-80,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(volSpec.y, volBC.u) annotation (Line(
          points={{-90,9},{-90,0},{-80,0},{-80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_N, volBC.u) annotation (Line(
          points={{-70,40},{-70,0},{-80,0},{-80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.N, u_N) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{-70,50},{-70,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 1st component of velocity
      connect(lin1BC.inert, inert) annotation (Line(
          points={{-40,-20},{-40,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin1Spec.y, lin1BC.u) annotation (Line(
          points={{-50,9},{-50,0},{-40,0},{-40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi_1, lin1BC.u) annotation (Line(
          points={{-30,40},{-30,0},{-40,0},{-40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi_1, u_Phi_1) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{-30,50},{-30,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 2nd component of velocity
      connect(lin2BC.inert, inert) annotation (Line(
          points={{6.10623e-16,-20},{6.10623e-16,-30},{5.55112e-16,-30},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin2Spec.y, lin2BC.u) annotation (Line(
          points={{-10,9},{-10,0},{0,0},{0,-12},{6.10623e-16,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi_2, lin2BC.u) annotation (Line(
          points={{10,40},{10,0},{6.10623e-16,0},{6.10623e-16,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi_2, u_Phi_2) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{10,50},{10,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 3rd component of velocity
      connect(lin3BC.inert, inert) annotation (Line(
          points={{40,-20},{40,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin3Spec.y, lin3BC.u) annotation (Line(
          points={{30,9},{30,0},{40,0},{40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi3, lin3BC.u) annotation (Line(
          points={{50,40},{50,0},{40,0},{40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi3, u_Phi3) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{50,50},{50,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(heatBC.inert, inert) annotation (Line(
          points={{80,-20},{80,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(heatSpec.y, heatBC.u) annotation (Line(
          points={{70,9},{70,0},{80,0},{80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Q, heatBC.u) annotation (Line(
          points={{90,40},{90,0},{80,0},{80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Q, u_Q) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{90,50},{90,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      annotation (
        defaultComponentName="phaseInertBC",
        Diagram(graphics),
        Icon(graphics));
    end Phase;

    package Volume "BCs for additivity of volume"
      extends Modelica.Icons.Package;

      model Volume "Prescribed volume"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Volume,
            redeclare Connectors.RealInput u(final unit="l3"));
      equation
        inert.V = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volBC");
      end Volume;

      model Pressure "Prescribed pressure"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Pressure,
            redeclare Connectors.RealInput u(final unit="m/(l.T2)"));
      equation
        inert.p = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volBC");
      end Pressure;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a material BC"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.mPhidot = zeros(n_vel) "No force";
          inert.Qdot = 0 "Adiabatic";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="volBC");
        end PartialBC;

        type BCType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of BCs";
      end BaseClasses;
    end Volume;

    package Momentum "BCs for linear momentum"
      extends Modelica.Icons.Package;
      model Velocity "Prescribed velocity"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            redeclare Connectors.RealInput u(final unit="l/T"));
      equation
        inert.phi[ax] = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="momBC");
      end Velocity;

      model Force "Prescribed force"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force,
            redeclare Connectors.RealInput u(final unit="l.m/T2"));
      equation
        inert.mPhidot[ax] = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="momBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for linear momentum"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;
          parameter Integer ax(
            final min=1,
            final max=n_vel) = 1
            "<html>Axis (between 1 and <i>n</i><sub>vel</sub>)</html>"
            annotation (HideResult=true);
          constant FCSys.BCs.InertAmagat.Momentum.BaseClasses.BCType bCType
            "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.V = 0 "No volume";
          for i in 1:n_vel loop
            if i <> ax then
              inert.mPhidot[i] = 0 "No force along the other axes";
            end if;
          end for;
          inert.Qdot = 0 "Adiabatic";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of BCs";
      end BaseClasses;
    end Momentum;

    package Heat "BCs for heat"
      extends Modelica.Icons.Package;

      model Temperature "Prescribed temperature"
        extends FCSys.BCs.InertAmagat.Heat.BaseClasses.PartialBC(final bCType=
              BaseClasses.BCType.Temperature, redeclare Connectors.RealInput u(
              final unit="l2.m/(N.T2)", displayUnit="K"));
      equation
        inert.T = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="heatBC");
      end Temperature;

      model HeatFlowRate "Prescribed heat flow rate"
        extends FCSys.BCs.InertAmagat.Heat.BaseClasses.PartialBC(final bCType=
              BaseClasses.BCType.HeatFlowRate, redeclare Connectors.RealInput u(
              final unit="l2.m/T3"));
      equation
        inert.Qdot = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="heatBC");
      end HeatFlowRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for heat"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;
          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        equation
          inert.V = 0 "No volume";
          inert.mPhidot = zeros(n_vel) "No force";
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="heatBC",
            Icon(graphics));
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Heat;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Integer n_vel(
          final min=1,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.InertAmagat inert(final n_vel=n_vel)
          "Connector for linear momentum and heat, with additivity of volume"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-50},
                  {10,-30}})));
        Connectors.RealInput u "Value of BC" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseInertBC");
      end PartialBC;
    end BaseClasses;
  end InertAmagat;

  package InertDalton
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>BC for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model</html>"

      extends FCSys.BaseClasses.Icons.BCs.Single;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>";

      // Additivity of pressure
      replaceable FCSys.BCs.InertDalton.Pressure.Volume pressBC(final n_vel=
            n_vel) constrainedby
        FCSys.BCs.InertDalton.Pressure.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Pressure", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-90,-26},{-70,-6}})));
      parameter Boolean internalPress=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Pressure",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant pressSpec(k(start=1*U.cm^3))
        if internalPress constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Pressure",
          __Dymola_descriptionLabel=true,
          enable=internalPress),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,20})));

      // 1st component of velocity
      replaceable FCSys.BCs.InertDalton.Momentum.Velocity lin1BC(final n_vel=
            n_vel) if n_vel > 0 constrainedby
        FCSys.BCs.InertDalton.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="1st component of velocity",
          enable=n_vel > 0,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-50,-26},{-30,-6}})));
      parameter Boolean internalLin1=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="1st component of velocity",
          __Dymola_descriptionLabel=true,
          enable=n_vel > 0,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin1Spec(k(start=0)) if
        internalLin1 and n_vel > 0 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal signal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="1st component of velocity",
          __Dymola_descriptionLabel=true,
          enable=n_vel > 0 and internalLin1),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,20})));

      // 2nd component of velocity
      replaceable FCSys.BCs.InertDalton.Momentum.Velocity lin2BC(final n_vel=
            n_vel) if n_vel > 1 constrainedby
        FCSys.BCs.InertDalton.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-26},{10,-6}})));
      parameter Boolean internalLin2=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin2Spec(k(start=0)) if
        internalLin2 and n_vel > 1 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="2nd component of velocity",
          enable=n_vel > 1 and internalLin2,
          __Dymola_descriptionLabel=true,
          enable=internalLin2),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,20})));

      // 3rd component of velocity
      replaceable FCSys.BCs.InertDalton.Momentum.Velocity lin3BC(final n_vel=
            n_vel) if n_vel > 2 constrainedby
        FCSys.BCs.InertDalton.Momentum.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{30,-26},{50,-6}})));
      parameter Boolean internalLin3=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant lin3Spec(k(start=0)) if
        internalLin3 and n_vel > 2 constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="3rd component of velocity",
          enable=n_vel > 2 and internalLin3,
          __Dymola_descriptionLabel=true,
          enable=internalLin3),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,20})));

      // Heat
      replaceable FCSys.BCs.InertDalton.Heat.Temperature heatBC(final n_vel=
            n_vel) constrainedby
        FCSys.BCs.InertDalton.Heat.BaseClasses.PartialBC "Type of condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Heat", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{70,-26},{90,-6}})));
      parameter Boolean internalHeat=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Heat",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant heatSpec(k(start=298.15*U.K))
        if internalHeat constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Heat",
          __Dymola_descriptionLabel=true,
          enable=internalHeat),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={70,20})));

      Connectors.InertDalton inert(final n_vel=n_vel)
        "Single-species connector for linear momentum and heat, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
      Connectors.RealInputBus u "Input bus for external signal sources"
        annotation (HideResult=not (internalPress or internalLin1 or
            internalLin2 or internalLin3 or internalHeat), Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,70}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

    protected
      Connectors.RealInputInternal u_N if not internalPress
        "External material signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,40})));
      Connectors.RealInputInternal u_Phi_1 if not internalLin1
        "External signal for 1st axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Phi_2 if not internalLin2 and n_vel > 1
        "External signal for 2nd axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Phi3 if not internalLin3 and n_vel > 2
        "External signal for 3rd axis of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,40})));
      Connectors.RealInputInternal u_Q if not internalHeat
        "External signal for heat" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,40}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,40})));

    equation
      // Material
      connect(pressBC.inert, inert) annotation (Line(
          points={{-80,-20},{-80,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(pressSpec.y, pressBC.u) annotation (Line(
          points={{-90,9},{-90,0},{-80,0},{-80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_N, pressBC.u) annotation (Line(
          points={{-70,40},{-70,0},{-80,0},{-80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.N, u_N) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{-70,50},{-70,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 1st axis
      connect(lin1BC.inert, inert) annotation (Line(
          points={{-40,-20},{-40,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin1Spec.y, lin1BC.u) annotation (Line(
          points={{-50,9},{-50,0},{-40,0},{-40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi_1, lin1BC.u) annotation (Line(
          points={{-30,40},{-30,0},{-40,0},{-40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi_1, u_Phi_1) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{-30,50},{-30,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 2nd axis
      connect(lin2BC.inert, inert) annotation (Line(
          points={{6.10623e-16,-20},{6.10623e-16,-30},{5.55112e-16,-30},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin2Spec.y, lin2BC.u) annotation (Line(
          points={{-10,9},{-10,0},{0,0},{0,-12},{6.10623e-16,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi_2, lin2BC.u) annotation (Line(
          points={{10,40},{10,0},{6.10623e-16,0},{6.10623e-16,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi_2, u_Phi_2) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{10,50},{10,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // 3rd axis
      connect(lin3BC.inert, inert) annotation (Line(
          points={{40,-20},{40,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(lin3Spec.y, lin3BC.u) annotation (Line(
          points={{30,9},{30,0},{40,0},{40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Phi3, lin3BC.u) annotation (Line(
          points={{50,40},{50,0},{40,0},{40,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Phi3, u_Phi3) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{50,50},{50,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(heatBC.inert, inert) annotation (Line(
          points={{80,-20},{80,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(heatSpec.y, heatBC.u) annotation (Line(
          points={{70,9},{70,0},{80,0},{80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_Q, heatBC.u) annotation (Line(
          points={{90,40},{90,0},{80,0},{80,-12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u.Q, u_Q) annotation (Line(
          points={{5.55112e-16,70},{0,70},{0,50},{90,50},{90,40}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      annotation (
        defaultComponentName="speciesInertBC",
        Diagram(graphics),
        Icon(graphics));
    end Species;

    package Pressure "BCs for additivity of pressure"
      extends Modelica.Icons.Package;

      model Volume "Prescribed volume"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Volume,
            redeclare Connectors.RealInput u(final unit="l3"));
      equation
        inert.V = u;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressBC");
      end Volume;

      model Pressure "Prescribed pressure"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Pressure,
            redeclare Connectors.RealInput u(final unit="m/(l.T2)"));
      equation
        inert.p = u;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressBC");
      end Pressure;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a material BC"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.mPhidot = zeros(n_vel) "No force";
          inert.Qdot = 0 "No heat flow";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="pressBC");
        end PartialBC;

        type BCType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of BCs";
      end BaseClasses;
    end Pressure;

    package Momentum "BCs for linear momentum"
      extends Modelica.Icons.Package;
      model Velocity "Prescribed velocity"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            redeclare Connectors.RealInput u(final unit="l/T"));
      equation
        inert.phi[ax] = u;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="momBC");
      end Velocity;

      model Force "Prescribed force"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force,
            redeclare Connectors.RealInput u(final unit="l.m/T2"));
      equation
        inert.mPhidot[ax] = u;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="momBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for linear momentum"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;
          parameter Integer ax(
            final min=1,
            final max=n_vel) = 1
            "<html>Axis (between 1 and <i>n</i><sub>vel</sub>)</html>"
            annotation (HideResult=true);
          constant FCSys.BCs.InertDalton.Momentum.BaseClasses.BCType bCType
            "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.p = 0 "No pressure";
          for i in 1:n_vel loop
            if i <> ax then
              inert.mPhidot[i] = 0 "No force along the other axes";
            end if;
          end for;
          inert.Qdot = 0 "Adiabatic";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of BCs";
      end BaseClasses;
    end Momentum;

    package Heat "BCs for heat"
      extends Modelica.Icons.Package;

      model Temperature "Prescribed temperature"
        extends FCSys.BCs.InertDalton.Heat.BaseClasses.PartialBC(final bCType=
              BaseClasses.BCType.Temperature, redeclare Connectors.RealInput u(
              final unit="l2.m/(N.T2)", displayUnit="K"));
      equation
        inert.T = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="heatBC");
      end Temperature;

      model HeatFlowRate "Prescribed heat flow rate"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.HeatFlowRate,
            redeclare Connectors.RealInput u(final unit="l2.m/T3"));
      equation
        inert.Qdot = u;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="heatBC");
      end HeatFlowRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for heat"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;
          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.p = 0 "No pressure";
          inert.mPhidot = zeros(n_vel) "No force";
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="heatBC",
            Icon(graphics));
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Heat;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Integer n_vel(
          final min=0,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);

        Connectors.InertDalton inert(final n_vel=n_vel)
          "Connector for linear momentum and heat, with additivity of pressure"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-50},
                  {10,-30}})));
        Connectors.RealInput u "Value of BC" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="speciesInertBC",
          Diagram(graphics));
      end PartialBC;
    end BaseClasses;
  end InertDalton;

  package Face
    "<html>BCs for <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> and <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"
    extends Modelica.Icons.Package;
    // **Update this package.

    model Subregion
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts by default</html>"
      import FCSys.BaseClasses.Axis;
      extends FCSys.BaseClasses.Icons.BCs.Single;

      parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

      FCSys.BCs.Face.Phases.Gas gas(final axis=axis) "Gas" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      FCSys.BCs.Face.Phases.Graphite graphite(final axis=axis) "Graphite"
        annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.BCs.Face.Phases.Ionomer ionomer(final axis=axis) "Ionomer"
        annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.FaceBus face
        "Connector for material, linear momentum, and heat of multiple species"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
      Connectors.RealInputBus u "Bus of inputs to specify conditions"
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

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregionFaceBC",
        Icon(graphics),
        Diagram(graphics));
    end Subregion;

    model SubregionFlow
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows by default</html>"
      extends FCSys.BCs.Face.Subregion(
        gas(
          H2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          N2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          O2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))),
        graphite(C(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)), 'e-'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))),
        ionomer(
          C19HF37O5S(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          'H+'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionFlow;

    model SubregionFlowTemp
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with current and force by default</html>"
      extends FCSys.BCs.Face.Subregion(
        gas(
          H2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC),
          N2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC),
          O2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC)),
        graphite(C(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC), 'e-'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC)),
        ionomer(
          C19HF37O5S(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC),
          'H+'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC)));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionFlowTemp;

    model Subregion0Current
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except zero current by default</html>"
      extends FCSys.BCs.Face.Subregion(
        gas(
          H2(redeclare replaceable Species.Material.Current materialBC),
          H2O(redeclare replaceable Species.Material.Current materialBC),
          N2(redeclare replaceable Species.Material.Current materialBC),
          O2(redeclare replaceable Species.Material.Current materialBC)),
        graphite(C(redeclare replaceable Species.Material.Current materialBC),
            'e-'(redeclare replaceable Species.Material.Current materialBC)),
        ionomer(
          C19HF37O5S(redeclare replaceable Species.Material.Current materialBC),

          H2O(redeclare replaceable Species.Material.Current materialBC),
          'H+'(redeclare replaceable Species.Material.Current materialBC)));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Current;

    model Subregion0Current0Power
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows except zero velocities by default</html>"
      extends FCSys.BCs.Face.Subregion(
        gas(
          H2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          N2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          O2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))),
        graphite(C(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))), 'e-'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))),
        ionomer(
          C19HF37O5S(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          'H+'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Current0Power;

    model Subregion0Power
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except zero power by default</html>"
      extends FCSys.BCs.Face.Subregion(
        gas(
          H2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0))),
          H2O(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0))),
          N2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0))),
          O2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0)))),
        graphite(C(redeclare replaceable Species.Heat.HeatFlowRate heatBC,
              heatSpec(k(start=0))), 'e-'(redeclare replaceable
              Species.Heat.HeatFlowRate heatBC, heatSpec(k(start=0)))),
        ionomer(
          C19HF37O5S(redeclare replaceable Species.Heat.HeatFlowRate heatBC,
              heatSpec(k(start=0))),
          H2O(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0))),
          'H+'(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Power;

    package Phases
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;
      model Phase "BC for a phase with all species conditionally included"

        extends FCSys.BCs.Face.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species C(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));
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
        FCSys.BCs.Face.Species.Species C19HF37O5S(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species 'e-' if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species H2 if inclH2 "Model" annotation (Dialog(
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
        FCSys.BCs.Face.Species.Species H2O if inclH2O "Model" annotation (
            Dialog(
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
        FCSys.BCs.Face.Species.Species 'H+' if 'inclH+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species N2 if inclN2 "Model" annotation (Dialog(
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
        FCSys.BCs.Face.Species.Species O2 if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // Note:  It would be helpful if Modelica allowed elements of expandable
        // connectors to be named by the contents of a string variable and the
        // name of an instance of a model was accessible through a string (like
        // %name is expanded to be the name of the instance of the model).  Then,
        // the connection equations that follow could be generic.

        // C
        connect(C.face.material, face.C.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.thermal, face.C.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, x.C.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, y.C.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, z.C.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, x.C.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, y.C.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, z.C.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.face.material, face.C19HF37O5S.material) annotation
          (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.thermal, face.C19HF37O5S.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, x.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, y.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, z.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, x.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, y.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, z.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.material, face.'e-'.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.thermal, face.'e-'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, x.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, y.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, z.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, x.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, y.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, z.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2
        connect(H2.face.material, face.H2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.thermal, face.H2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, x.H2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, y.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, z.H2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, x.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, y.H2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, z.H2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect(H2O.face.material, face.H2O.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, x.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, y.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, z.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, x.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, y.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, z.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect('H+'.face.material, face.'H+'.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, x.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, y.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, z.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, x.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, y.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, z.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face.material, face.N2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.thermal, face.N2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, x.N2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, y.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, z.N2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, x.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, y.N2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, z.N2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect(O2.face.material, face.O2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.thermal, face.O2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, x.O2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, y.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, z.O2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, x.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, y.O2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, z.O2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Phase;

      model Gas "BC for gas"

        extends FCSys.BCs.Face.Phases.BaseClasses.NullPhase;

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
        FCSys.BCs.Face.Species.Species H2 if inclH2 "Model" annotation (Dialog(
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
        FCSys.BCs.Face.Species.Species H2O if inclH2O "Model" annotation (
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
        FCSys.BCs.Face.Species.Species N2 if inclN2 "Model" annotation (Dialog(
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
        FCSys.BCs.Face.Species.Species O2 if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.face.material, face.H2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.thermal, face.H2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, x.H2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, y.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical1, z.H2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, x.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, y.H2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.mechanical2, z.H2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect(H2O.face.material, face.H2O.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, x.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, y.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical1, z.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, x.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, y.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.mechanical2, z.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect(N2.face.material, face.N2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.thermal, face.N2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, x.N2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, y.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical1, z.N2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, x.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, y.N2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.mechanical2, z.N2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
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
        connect(O2.face.material, face.O2.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.thermal, face.O2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, x.O2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, y.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical1, z.O2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, x.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, y.O2.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.mechanical2, z.O2.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Gas;

      model Graphite "BC for graphite"

        extends FCSys.BCs.Face.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species C(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.Face.Species.Species 'e-' if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C
        connect(C.face.material, face.C.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.thermal, face.C.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, x.C.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, y.C.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical1, z.C.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, x.C.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, y.C.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.mechanical2, z.C.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.material, face.'e-'.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.thermal, face.'e-'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, x.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, y.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical1, z.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, x.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, y.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.mechanical2, z.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Graphite;

      model Ionomer "BC for ionomer"

        extends FCSys.BCs.Face.Phases.BaseClasses.NullPhase;

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
        FCSys.BCs.Face.Species.Species C19HF37O5S(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
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
        FCSys.BCs.Face.Species.Species H2O if inclH2O "Model" annotation (
            Dialog(
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
        FCSys.BCs.Face.Species.Species 'H+' if 'inclH+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.face.material, face.C19HF37O5S.material) annotation
          (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.thermal, face.C19HF37O5S.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, x.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, y.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical1, z.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, x.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, y.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.mechanical2, z.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face.material, face.'H+'.material) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, x.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, y.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical1, z.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, x.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{-20,-10},{-20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, y.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.mechanical2, z.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,-4},{0,-10},{20,-10},{20,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceBC",
          Diagram(graphics));
      end Ionomer;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty BC for a phase (no species)"
          import FCSys.BaseClasses.Axis;
          extends FCSys.BaseClasses.Icons.BCs.Single;

          parameter FCSys.BaseClasses.Axis axis=Axis.x
            "Axis normal to the face";

          FCSys.Connectors.FaceBus face
            "Multi-species connector for material, linear momentum, and heat"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,40})));

        protected
          Connectors.FaceBusInternal x if axis == Axis.x
            "Internal connector enabled if x axis" annotation (Placement(
                transformation(extent={{-30,-30},{-10,-10}})));
          Connectors.FaceBusInternal y if axis == Axis.y
            "Internal connector enabled if y axis"
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
          Connectors.FaceBusInternal z if axis == Axis.z
            "Internal connector enabled if z axis"
            annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
        equation

          connect(x, face) annotation (Line(
              points={{-20,-20},{-20,-30},{5.55112e-16,-30},{5.55112e-16,-40}},

              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(y, face) annotation (Line(
              points={{5.55112e-16,-20},{5.55112e-16,-30},{5.55112e-16,-40},{
                  5.55112e-16,-40}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(z, face) annotation (Line(
              points={{20,-20},{20,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC",
            Diagram(graphics));
        end NullPhase;
      end BaseClasses;
    end Phases;

    package Species
      "<html>BCs for a single <a href=\"modelica://FCSys.Connectors.BaseClasses.BaseClasses.PartialFace\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends Modelica.Icons.Package;

      model Species
        "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter ThermoOpt thermoOpt=ThermoOpt.OpenDiabatic
          "Options for material and thermal transport";

        // Material
        final parameter Boolean open=thermoOpt == ThermoOpt.OpenDiabatic "Open";
        // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
        // option, e.g.,
        //     enable=thermoOpt=ThermoOpt.OpenDiabatic.
        // Therefore, the values of the enumerations are specified numerically for
        // this initial condition and others below for material and heat.
        replaceable FCSys.BCs.Face.Species.Material.Pressure materialBC if open
          constrainedby FCSys.BCs.Face.Species.Material.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-70,-26},{-50,-6}})));

        parameter Boolean internalMaterial=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant materialSpec(k(start=0))
          if internalMaterial and open constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            enable=thermoOpt == 3 and internalMaterial,
            __Dymola_descriptionLabel=true,
            enable=thermoOpt == 3 and internalMaterial),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-70,20})));

        // 1st transverse linear momentum
        parameter Boolean slip1=false
          "<html>Viscous (1<sup>st</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="1st transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Momentum.Velocity lin1BC if slip1 constrainedby
          FCSys.BCs.Face.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-30,-26},{-10,-6}})));
        parameter Boolean internalLin1=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Modelica.Blocks.Sources.Constant lin1Spec(k(start=0)) if
          internalLin1 and slip1 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=internalLin1 and slip1,
            __Dymola_descriptionLabel=true,
            enable=internalLin1 and slip1),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,20})));

        // 2nd transverse linear momentum
        parameter Boolean slip2=false
          "<html>Viscous (2<sup>nd</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="2nd transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Momentum.Velocity lin2BC if slip2 constrainedby
          FCSys.BCs.Face.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{10,-26},{30,-6}})));

        parameter Boolean internalLin2=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant lin2Spec(k(start=0)) if
          internalLin2 and slip2 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=internalLin2 and slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,20})));

        // Heat
        final parameter Boolean diabatic=thermoOpt == ThermoOpt.ClosedDiabatic
             or thermoOpt == ThermoOpt.OpenDiabatic "Diabatic (heat included)";
        replaceable Heat.Temperature heatBC if diabatic constrainedby
          FCSys.BCs.Face.Species.Heat.BaseClasses.PartialBC "Type of condition"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{50,-26},{70,-6}})));

        parameter Boolean internalHeat=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Modelica.Blocks.Sources.Constant heatSpec(k(start=298.15*U.K))
          if internalHeat and diabatic constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            __Dymola_descriptionLabel=true,
            enable=(thermoOpt == 2 or thermoOpt == 3) and internalHeat,
            enable=internalHeat and diabatic),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={50,20})));

        FCSys.Connectors.FaceGeneric face(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Single-species connector for material, linear momentum, and heat"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-50},{10,-30}})));
        FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalHeat), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,70}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        Connectors.RealInputInternal u_N if not internalMaterial and open
          "External material signal" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-50,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-60,40})));
        Connectors.RealInputInternal u_Phi_1 if not internalLin1 and slip1
          "External signal for 1st transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Phi_2 if not internalLin2 and slip2
          "External signal for 2nd transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Q if not internalHeat and diabatic
          "External signal for heat" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,40})));

      public
        Connectors.RealInput u1 "**test" annotation (Dialog, Placement(
              transformation(extent={{-88,78},{-68,98}})));
      equation
        // Material
        connect(materialBC.material, face.material) annotation (Line(
            points={{-60,-20},{-60,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(materialSpec.y, materialBC.u) annotation (Line(
            points={{-70,9},{-70,0},{-60,0},{-60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_N, materialBC.u) annotation (Line(
            points={{-50,40},{-50,0},{-60,0},{-60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.N, u_N) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-50,50},{-50,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 1st transverse linear momentum
        connect(lin1BC.momentum, face.mechanical1) annotation (Line(
            points={{-20,-20},{-20,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin1Spec.y, lin1BC.u) annotation (Line(
            points={{-30,9},{-30,0},{-20,0},{-20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_1, lin1BC.u) annotation (Line(
            points={{-10,40},{-10,0},{-20,0},{-20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_1, u_Phi_1) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-10,50},{-10,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 2nd transverse linear momentum
        connect(lin2BC.momentum, face.mechanical2) annotation (Line(
            points={{20,-20},{20,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin2Spec.y, lin2BC.u) annotation (Line(
            points={{10,9},{10,0},{20,0},{20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_2, lin2BC.u) annotation (Line(
            points={{30,40},{30,0},{20,0},{20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_2, u_Phi_2) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{30,50},{30,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Heat
        connect(heatBC.thermal, face.thermal) annotation (Line(
            points={{60,-20},{60,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(heatSpec.y, heatBC.u) annotation (Line(
            points={{50,9},{50,0},{60,0},{60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Q, heatBC.u) annotation (Line(
            points={{70,40},{70,0},{60,0},{60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Q, u_Q) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{70,50},{70,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        annotation (
          defaultComponentName="speciesFaceBC",
          Diagram(graphics),
          Icon(graphics));
      end Species;

      package Material "BCs for material"
        extends Modelica.Icons.Package;

        model Pressure "Prescribed pressure"

          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Pressure,
              redeclare Connectors.RealInput u(final unit="m/(l.T2)"));
          // **Press->EC pot
        equation
          material.mu = u;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="materialBC",
            Documentation(info="<html><p>Pressure can be
  calculated from other properties (e.g., temperature, specific volume, specific enthalpy, specific heat, or Gibbs potential)
  using functions in the
  <a href=\"modelica://FCSys.Connectors.Characteristic\">Characteristics</a> package.
  </p></html>"));
        end Pressure;

        model Current "Prescribed current"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Current,
              redeclare Connectors.RealInput u(final unit="N/T"));

        equation
          material.Ndot = u;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="materialBC");
        end Current;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a material BC"
            extends FCSys.BaseClasses.Icons.BCs.Single;
            constant BCType bCType "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.MaterialTransport material
              "Material connector for the face" annotation (Placement(
                  transformation(extent={{-10,-50},{10,-30}})));
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="materialBC",
              Diagram(graphics));
          end PartialBC;

          type BCType = enumeration(
              Pressure "Pressure",
              Current "Current") "Types of BCs";
        end BaseClasses;
      end Material;

      package Momentum "BCs for linear momentum"
        extends Modelica.Icons.Package;
        model Velocity "Prescribed velocity"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
              redeclare Connectors.RealInput u(final unit="l/T"));
        equation
          momentum.phi = u;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end Velocity;

        model Force "Prescribed force"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force,
              redeclare Connectors.RealInput u(final unit="l.m/T2"));
        equation
          momentum.mPhidot = u;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end Force;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a BC for linear momentum"
            extends FCSys.BaseClasses.Icons.BCs.Single;
            constant FCSys.BCs.Face.Species.Momentum.BaseClasses.BCType bCType
              "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.MechanicalTransport momentum
              "Linear momentum connector for the face" annotation (Placement(
                  transformation(extent={{-10,-50},{10,-30}})));
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="momBC");
          end PartialBC;

          type BCType = enumeration(
              Velocity "Velocity",
              Force "Force") "Types of BCs";
        end BaseClasses;
      end Momentum;

      package Heat "BCs for heat"
        extends Modelica.Icons.Package;

        model Temperature "Prescribed temperature"
          extends FCSys.BCs.Face.Species.Heat.BaseClasses.PartialBC(final
              bCType=BaseClasses.BCType.Temperature, redeclare
              Connectors.RealInput u(final unit="l2.m/(N.T2)", displayUnit="K"));
        equation
          heat.T = u;
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="heatBC");
        end Temperature;

        model HeatFlowRate "Prescribed heat flow rate"
          extends FCSys.BCs.Face.Species.Heat.BaseClasses.PartialBC(final
              bCType=BaseClasses.BCType.HeatFlowRate, redeclare
              Connectors.RealInput u(final unit="l2.m/T3"));

        equation
          heat.Qdot = u;
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="heatBC",
            Diagram(graphics));
        end HeatFlowRate;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a BC for heat"
            extends FCSys.BaseClasses.Icons.BCs.Single;
            constant BCType bCType "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.Thermal thermal "Thermal connector for the face"
              annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          equation

            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="heatBC",
              Diagram(graphics));
          end PartialBC;

          type BCType = enumeration(
              Temperature "Temperature",
              HeatFlowRate "Heat flow rate") "Types of BCs";
        end BaseClasses;
      end Heat;

      model SpeciesX
        "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        // **Create models for Generic, Y, and Z faces.
        // Material
        final parameter Boolean open=thermoOpt == ThermoOpt.OpenDiabatic "Open";
        // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
        // option, e.g.,
        //     enable=thermoOpt=ThermoOpt.OpenDiabatic.
        // Therefore, the values of the enumerations are specified numerically for
        // this initial condition and others below for material and heat.
        replaceable FCSys.BCs.Face.Species.Material.Pressure materialBC if open
          constrainedby FCSys.BCs.Face.Species.Material.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-70,-26},{-50,-6}})));

        parameter Boolean internalMaterial=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant materialSpec(k(start=0))
          if internalMaterial and open constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            enable=thermoOpt == 3 and internalMaterial,
            __Dymola_descriptionLabel=true,
            enable=thermoOpt == 3 and internalMaterial),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-70,20})));

        // 1st transverse linear momentum
        parameter Boolean slip1=false
          "<html>Viscous (1<sup>st</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="1st transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable FCSys.BCs.Face.Species.Momentum.Velocity lin1BC if slip1
          constrainedby FCSys.BCs.Face.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-30,-26},{-10,-6}})));
        parameter Boolean internalLin1=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Modelica.Blocks.Sources.Constant lin1Spec(k(start=0)) if
          internalLin1 and slip1 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=internalLin1 and slip1,
            __Dymola_descriptionLabel=true,
            enable=internalLin1 and slip1),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,20})));

        // 2nd transverse linear momentum
        parameter Boolean slip2=false
          "<html>Viscous (2<sup>nd</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="2nd transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable FCSys.BCs.Face.Species.Momentum.Velocity lin2BC if slip2
          constrainedby FCSys.BCs.Face.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{10,-26},{30,-6}})));

        parameter Boolean internalLin2=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant lin2Spec(k(start=0)) if
          internalLin2 and slip2 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=internalLin2 and slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,20})));

        // Heat
        final parameter Boolean diabatic=thermoOpt == ThermoOpt.ClosedDiabatic
             or thermoOpt == ThermoOpt.OpenDiabatic "Diabatic (heat included)";
        replaceable FCSys.BCs.Face.Species.Heat.Temperature heatBC if diabatic
          constrainedby FCSys.BCs.Face.Species.Heat.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{50,-26},{70,-6}})));

        parameter Boolean internalHeat=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Modelica.Blocks.Sources.Constant heatSpec(k(start=298.15*U.K))
          if internalHeat and diabatic constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            __Dymola_descriptionLabel=true,
            enable=(thermoOpt == 2 or thermoOpt == 3) and internalHeat,
            enable=internalHeat and diabatic),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={50,20})));

        Connectors.FaceX face(final thermoOpt=thermoOpt)
          "Single-species connector for material, linear momentum, and heat"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-50},{10,-30}})));
        FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalHeat), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,70}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        Connectors.RealInputInternal u_N if not internalMaterial and open
          "External material signal" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-50,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-60,40})));
        Connectors.RealInputInternal u_Phi_1 if not internalLin1 and slip1
          "External signal for 1st transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Phi_2 if not internalLin2 and slip2
          "External signal for 2nd transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Q if not internalHeat and diabatic
          "External signal for heat" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,40})));

      equation
        // Material
        connect(materialBC.material, face.material) annotation (Line(
            points={{-60,-20},{-60,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(materialSpec.y, materialBC.u) annotation (Line(
            points={{-70,9},{-70,0},{-60,0},{-60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_N, materialBC.u) annotation (Line(
            points={{-50,40},{-50,0},{-60,0},{-60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.N, u_N) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-50,50},{-50,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 1st transverse linear momentum
        connect(lin1Spec.y, lin1BC.u) annotation (Line(
            points={{-30,9},{-30,0},{-20,0},{-20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_1, lin1BC.u) annotation (Line(
            points={{-10,40},{-10,0},{-20,0},{-20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_1, u_Phi_1) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-10,50},{-10,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 2nd transverse linear momentum
        connect(lin2Spec.y, lin2BC.u) annotation (Line(
            points={{10,9},{10,0},{20,0},{20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_2, lin2BC.u) annotation (Line(
            points={{30,40},{30,0},{20,0},{20,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_2, u_Phi_2) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{30,50},{30,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Heat
        connect(heatBC.thermal, face.thermal) annotation (Line(
            points={{60,-20},{60,-30},{5.55112e-16,-30},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(heatSpec.y, heatBC.u) annotation (Line(
            points={{50,9},{50,0},{60,0},{60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Q, heatBC.u) annotation (Line(
            points={{70,40},{70,0},{60,0},{60,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Q, u_Q) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{70,50},{70,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        annotation (
          defaultComponentName="speciesFaceBC",
          Diagram(graphics),
          Icon(graphics));
      end SpeciesX;
    end Species;

    model Region

      Connectors.FaceBus face[1, 1]
        annotation (Placement(transformation(extent={{-32,-8},{-12,12}})));
      replaceable Subregion subregionFaceBC[1, 1](each graphite('incle-'=true,
            'e-'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC,
              redeclare Modelica.Blocks.Sources.Ramp materialSpec(duration=50,
                height=1*U.V))))
        annotation (Placement(transformation(extent={{32,-4},{52,16}})));
    equation
      connect(subregionFaceBC.face, face) annotation (Line(
          points={{42,2},{26,2},{26,2},{10,2},{10,2},{-22,2}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (Diagram(graphics));
    end Region;
    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the boundary conditions must be as well.  A
<a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a>
connector (<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
<a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>,
<a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>,
or <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a>)
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
  end Face;

  package FaceDifferential
    "<html>BCs for pairs of <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> or <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"
    extends Modelica.Icons.Package;
    // **Update this package.
    model Subregion
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts by default</html>"
      extends FCSys.BaseClasses.Icons.BCs.Double;

      replaceable FCSys.BCs.FaceDifferential.Phases.Gas gas(final axis=axis)
        "Gas" annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      replaceable FCSys.BCs.FaceDifferential.Phases.Graphite graphite(final
          axis=axis) "Graphite" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      replaceable FCSys.BCs.FaceDifferential.Phases.Ionomer ionomer(final axis=
            axis) "Ionomer" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      FCSys.Connectors.FaceBus negative annotation (Placement(transformation(
              extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      Connectors.RealInputBus u annotation (Placement(transformation(
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
        Icon(graphics),
        Diagram(graphics));
    end Subregion;

    model SubregionFlow
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows by default</html>"
      extends FCSys.BCs.FaceDifferential.Subregion(
        gas(
          H2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          N2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          O2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))),
        graphite(C(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)), 'e-'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))),
        ionomer(
          C19HF37O5S(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0)),
          'H+'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Momentum.Force lin1BC,
            redeclare replaceable Species.Momentum.Force lin2BC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k=0))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionFlow;

    model Subregion0Current
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except zero current by default</html>"
      extends FCSys.BCs.FaceDifferential.Subregion(
        gas(
          H2(redeclare replaceable Species.Material.Current materialBC),
          H2O(redeclare replaceable Species.Material.Current materialBC),
          N2(redeclare replaceable Species.Material.Current materialBC),
          O2(redeclare replaceable Species.Material.Current materialBC)),
        graphite(C(redeclare replaceable Species.Material.Current materialBC),
            'e-'(redeclare replaceable Species.Material.Current materialBC)),
        ionomer(
          C19HF37O5S(redeclare replaceable Species.Material.Current materialBC),

          H2O(redeclare replaceable Species.Material.Current materialBC),
          'H+'(redeclare replaceable Species.Material.Current materialBC)));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Current;

    model Subregion0Current0Power
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows except zero velocities by default</html>"
      extends FCSys.BCs.FaceDifferential.Subregion(
        gas(
          H2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          N2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          O2(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))),
        graphite(C(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))), 'e-'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))),
        ionomer(
          C19HF37O5S(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          H2O(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0))),
          'H+'(
            redeclare replaceable Species.Material.Current materialBC,
            redeclare replaceable Species.Heat.HeatFlowRate heatBC,
            heatSpec(k(start=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Current0Power;

    model Subregion0Power
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except zero power by default</html>"
      extends FCSys.BCs.FaceDifferential.Subregion(
        gas(
          H2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0))),
          H2O(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0))),
          N2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0))),
          O2(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(k(
                  start=0)))),
        graphite(C(redeclare replaceable Species.Heat.HeatFlowRate heatBC,
              heatSpec(k(start=0))), 'e-'(redeclare replaceable
              Species.Heat.HeatFlowRate heatBC, heatSpec(k(start=0)))),
        ionomer(
          C19HF37O5S(redeclare replaceable Species.Heat.HeatFlowRate heatBC,
              heatSpec(k(start=0))),
          H2O(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0))),
          'H+'(redeclare replaceable Species.Heat.HeatFlowRate heatBC, heatSpec(
                k(start=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end Subregion0Power;

    package Phases
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;
      model Phase "BC for a phase with all species conditionally included"

        extends FCSys.BCs.FaceDifferential.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=true "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species C(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));
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
        FCSys.BCs.FaceDifferential.Species.Species C19HF37O5S(thermoOpt=
              ThermoOpt.ClosedDiabatic) if inclC19HF37O5S "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species 'e-' if 'incle-' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species H2 if inclH2 "Model"
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
        FCSys.BCs.FaceDifferential.Species.Species H2O if inclH2O "Model"
          annotation (Dialog(
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
        FCSys.BCs.FaceDifferential.Species.Species 'H+' if 'inclH+' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species N2 if inclN2 "Model"
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
        FCSys.BCs.FaceDifferential.Species.Species O2 if inclO2 "Model"
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

        // C
        connect(C.negative.material, negative.C.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.thermal, negative.C.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, xNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, yNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, zNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, xNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, yNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, zNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.material, positive.C.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.thermal, positive.C.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, xPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, yPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, zPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, xPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, yPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, zPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.negative.material, negative.C19HF37O5S.material)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.thermal, negative.C19HF37O5S.thermal) annotation
          (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, xNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, yNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, zNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, xNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, yNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, zNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.material, positive.C19HF37O5S.material)
          annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.thermal, positive.C19HF37O5S.thermal) annotation
          (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, xPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, yPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, zPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, xPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, yPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, zPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative.material, negative.'e-'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.thermal, negative.'e-'.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, xNegative.'e-'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, yNegative.'e-'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, zNegative.'e-'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, xNegative.'e-'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, yNegative.'e-'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, zNegative.'e-'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.material, positive.'e-'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.thermal, positive.'e-'.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, xPositive.'e-'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, yPositive.'e-'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, zPositive.'e-'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, xPositive.'e-'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, yPositive.'e-'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, zPositive.'e-'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2
        connect(H2.negative.material, negative.H2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.thermal, negative.H2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, xNegative.H2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, yNegative.H2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, zNegative.H2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, xNegative.H2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, yNegative.H2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, zNegative.H2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.material, positive.H2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.thermal, positive.H2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, xPositive.H2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, yPositive.H2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, zPositive.H2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, xPositive.H2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, yPositive.H2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, zPositive.H2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect(H2O.negative.material, negative.H2O.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.thermal, negative.H2O.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, xNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, yNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, zNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, xNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, yNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, zNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.material, positive.H2O.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.thermal, positive.H2O.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, xPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, yPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, zPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, xPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, yPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, zPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect('H+'.negative.material, negative.'H+'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.thermal, negative.'H+'.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, xNegative.'H+'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, yNegative.'H+'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, zNegative.'H+'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, xNegative.'H+'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, yNegative.'H+'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, zNegative.'H+'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.material, positive.'H+'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.thermal, positive.'H+'.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, xPositive.'H+'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, yPositive.'H+'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, zPositive.'H+'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, xPositive.'H+'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, yPositive.'H+'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, zPositive.'H+'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.negative.material, negative.N2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.thermal, negative.N2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, xNegative.N2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, yNegative.N2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, zNegative.N2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, xNegative.N2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, yNegative.N2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, zNegative.N2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.material, positive.N2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.thermal, positive.N2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, xPositive.N2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, yPositive.N2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, zPositive.N2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, xPositive.N2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, yPositive.N2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, zPositive.N2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect(O2.negative.material, negative.O2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.thermal, negative.O2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, xNegative.O2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, yNegative.O2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, zNegative.O2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, xNegative.O2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, yNegative.O2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, zNegative.O2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.material, positive.O2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.thermal, positive.O2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, xPositive.O2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, yPositive.O2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, zPositive.O2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, xPositive.O2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, yPositive.O2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, zPositive.O2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Phase;

      model Gas "BC for gas"

        extends FCSys.BCs.FaceDifferential.Phases.BaseClasses.NullPhase;

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
        FCSys.BCs.FaceDifferential.Species.Species H2 if inclH2 "Model"
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
        FCSys.BCs.FaceDifferential.Species.Species H2O if inclH2O "Model"
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
        FCSys.BCs.FaceDifferential.Species.Species N2 if inclN2 "Model"
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
        FCSys.BCs.FaceDifferential.Species.Species O2 if inclO2 "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.negative.material, negative.H2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.thermal, negative.H2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, xNegative.H2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, yNegative.H2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical1, zNegative.H2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, xNegative.H2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, yNegative.H2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.mechanical2, zNegative.H2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.material, positive.H2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.thermal, positive.H2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, xPositive.H2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, yPositive.H2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical1, zPositive.H2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, xPositive.H2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, yPositive.H2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.mechanical2, zPositive.H2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect(H2O.negative.material, negative.H2O.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.thermal, negative.H2O.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, xNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, yNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, zNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, xNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, yNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, zNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.material, positive.H2O.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.thermal, positive.H2O.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, xPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, yPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, zPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, xPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, yPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, zPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect(N2.negative.material, negative.N2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.thermal, negative.N2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, xNegative.N2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, yNegative.N2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical1, zNegative.N2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, xNegative.N2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, yNegative.N2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.mechanical2, zNegative.N2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.material, positive.N2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.thermal, positive.N2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, xPositive.N2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, yPositive.N2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical1, zPositive.N2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, xPositive.N2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, yPositive.N2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.mechanical2, zPositive.N2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect(O2.negative.material, negative.O2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.thermal, negative.O2.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, xNegative.O2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, yNegative.O2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical1, zNegative.O2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-10,0},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, xNegative.O2.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, yNegative.O2.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.mechanical2, zNegative.O2.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.material, positive.O2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.thermal, positive.O2.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, xPositive.O2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, yPositive.O2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical1, zPositive.O2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, xPositive.O2.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, yPositive.O2.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.mechanical2, zPositive.O2.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Gas;

      model Graphite "BC for graphite"

        extends FCSys.BCs.FaceDifferential.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species C(thermoOpt=ThermoOpt.ClosedDiabatic)
          if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.BCs.FaceDifferential.Species.Species 'e-' if 'incle-' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C
        connect(C.negative.material, negative.C.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.thermal, negative.C.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, xNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, yNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical1, zNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, xNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, yNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.mechanical2, zNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.material, positive.C.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.thermal, positive.C.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, xPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, yPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical1, zPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, xPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, yPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.mechanical2, zPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C, C.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative.material, negative.'e-'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.thermal, negative.'e-'.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, xNegative.'e-'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, yNegative.'e-'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical1, zNegative.'e-'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, xNegative.'e-'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, yNegative.'e-'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.mechanical2, zNegative.'e-'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.material, positive.'e-'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.thermal, positive.'e-'.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, xPositive.'e-'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, yPositive.'e-'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical1, zPositive.'e-'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, xPositive.'e-'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, yPositive.'e-'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.mechanical2, zPositive.'e-'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC");
      end Graphite;

      model Ionomer "BC for ionomer"

        extends FCSys.BCs.FaceDifferential.Phases.BaseClasses.NullPhase;

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
        FCSys.BCs.FaceDifferential.Species.Species C19HF37O5S(thermoOpt=
              ThermoOpt.ClosedDiabatic) if inclC19HF37O5S "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
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
        FCSys.BCs.FaceDifferential.Species.Species H2O if inclH2O "Model"
          annotation (Dialog(
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
        FCSys.BCs.FaceDifferential.Species.Species 'H+' if 'inclH+' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.negative.material, negative.C19HF37O5S.material)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.thermal, negative.C19HF37O5S.thermal) annotation
          (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, xNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, yNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical1, zNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, xNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, yNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.mechanical2, zNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.material, positive.C19HF37O5S.material)
          annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.thermal, positive.C19HF37O5S.thermal) annotation
          (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, xPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, yPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical1, zPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, xPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, yPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.mechanical2, zPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.C19HF37O5S, C19HF37O5S.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative.material, negative.H2O.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.thermal, negative.H2O.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, xNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, yNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical1, zNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, xNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, yNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.mechanical2, zNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.material, positive.H2O.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.thermal, positive.H2O.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, xPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, yPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical1, zPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, xPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, yPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.mechanical2, zPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
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
        connect('H+'.negative.material, negative.'H+'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.thermal, negative.'H+'.thermal) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, xNegative.'H+'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, yNegative.'H+'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical1, zNegative.'H+'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, xNegative.'H+'.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, yNegative.'H+'.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.mechanical2, zNegative.'H+'.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.material, positive.'H+'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.thermal, positive.'H+'.thermal) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, xPositive.'H+'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, yPositive.'H+'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical1, zPositive.'H+'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, xPositive.'H+'.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, yPositive.'H+'.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.mechanical2, zPositive.'H+'.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceBC",
          Diagram(graphics));
      end Ionomer;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty BC for a phase (no species)"
          import FCSys.BaseClasses.Axis;
          extends FCSys.BaseClasses.Icons.BCs.Double;

          parameter FCSys.BaseClasses.Axis axis=Axis.x
            "Axis normal to the face";

          FCSys.Connectors.FaceBus negative
            "Multi-species connector for material, linear momentum, and heat"
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.FaceBus positive
            "Multi-species connectors for of material, linear momentum, and heat"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));

          Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,40})));

        protected
          Connectors.FaceBusInternal xNegative if axis == Axis.x
            "Internal connector for the negative face enabled if x axis"
            annotation (Placement(transformation(extent={{-70,10},{-50,30}})));
          Connectors.FaceBusInternal yNegative if axis == Axis.y
            "Internal connector for the negative face enabled if y axis"
            annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
          Connectors.FaceBusInternal zNegative if axis == Axis.z
            "Internal connector for the negative face enabled if z axis"
            annotation (Placement(transformation(extent={{-70,-30},{-50,-10}})));
          Connectors.FaceBusInternal xPositive if axis == Axis.x
            "Internal connector for the positive face enabled if x axis"
            annotation (Placement(transformation(extent={{50,10},{70,30}})));
          Connectors.FaceBusInternal yPositive if axis == Axis.y
            "Internal connector for the positive face enabled if y axis"
            annotation (Placement(transformation(extent={{50,-10},{70,10}})));
          Connectors.FaceBusInternal zPositive if axis == Axis.z
            "Internal connector for the positive face enabled if z axis"
            annotation (Placement(transformation(extent={{50,-30},{70,-10}})));
        equation

          connect(xNegative, negative) annotation (Line(
              points={{-60,20},{-80,20},{-80,0},{-100,0},{-100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(yNegative, negative) annotation (Line(
              points={{-60,5.55112e-16},{-100,0},{-100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(zNegative, negative) annotation (Line(
              points={{-60,-20},{-80,-20},{-80,0},{-100,0},{-100,5.55112e-16}},

              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(xPositive, positive) annotation (Line(
              points={{60,20},{80,20},{80,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          connect(yPositive, positive) annotation (Line(
              points={{60,5.55112e-16},{80,-4.87687e-22},{80,5.55112e-16},{100,
                  5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(zPositive, positive) annotation (Line(
              points={{60,-20},{80,-20},{80,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceBC",
            Diagram(graphics));
        end NullPhase;
      end BaseClasses;
    end Phases;

    package Species
      "<html>BCs for a Single <a href=\"modelica://FCSys.Connectors.BaseClasses.BaseClasses.PartialFace\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends Modelica.Icons.Package;

      model Species
        "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        extends FCSys.BaseClasses.Icons.BCs.Double;

        parameter ThermoOpt thermoOpt=ThermoOpt.OpenDiabatic
          "Options for material and thermal transport";

        // Material
        final parameter Boolean open=thermoOpt == ThermoOpt.OpenDiabatic "Open";
        replaceable
          FCSys.BCs.FaceDifferential.Species.Material.PotentialElectrochemical
          materialBC if open constrainedby
          FCSys.BCs.FaceDifferential.Species.Material.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-70,-30},{-50,-10}})));

        parameter Boolean internalMaterial=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Material",
            enable=thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant materialSpec(k(start=0))
          if internalMaterial and open constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Material",
            __Dymola_descriptionLabel=true,
            enable=thermoOpt == 3 and internalMaterial),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-70,20})));

        // 1st transverse linear momentum
        parameter Boolean slip1=false
          "<html>Viscous (1<sup>st</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="1st transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable FCSys.BCs.FaceDifferential.Species.Momentum.Velocity lin1BC
          if slip1 constrainedby
          FCSys.BCs.FaceDifferential.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-30,-40},{-10,-20}})));

        parameter Boolean internalLin1=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="1st transverse momentum",
            enable=slip1,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable Modelica.Blocks.Sources.Constant lin1Spec(k(start=0)) if
          internalLin1 and slip1 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="1st transverse momentum",
            enable=internalLin1 and slip1,
            __Dymola_descriptionLabel=true,
            enable=internalLin1 and slip1),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,20})));

        // 2nd transverse linear momentum
        parameter Boolean slip2=false
          "<html>Viscous (2<sup>nd</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="2nd transverse momentum",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));

        replaceable FCSys.BCs.FaceDifferential.Species.Momentum.Velocity lin2BC
          if slip2 constrainedby
          FCSys.BCs.FaceDifferential.Species.Momentum.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{10,-50},{30,-30}})));

        parameter Boolean internalLin2=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="2nd transverse momentum",
            enable=slip2,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant lin2Spec(k(start=0)) if
          internalLin2 and slip2 constrainedby Modelica.Blocks.Interfaces.SO
          "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="2nd transverse momentum",
            enable=internalLin2 and slip2,
            __Dymola_descriptionLabel=true),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,20})));

        // Heat
        final parameter Boolean diabatic=thermoOpt == ThermoOpt.ClosedDiabatic
             or thermoOpt == ThermoOpt.OpenDiabatic "Diabatic (heat included)";
        replaceable FCSys.BCs.FaceDifferential.Species.Heat.Temperature heatBC
          if diabatic constrainedby
          FCSys.BCs.FaceDifferential.Species.Heat.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{50,-60},{70,-40}})));
        parameter Boolean internalHeat=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Heat",
            enable=thermoOpt == 2 or thermoOpt == 3,
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Modelica.Blocks.Sources.Constant heatSpec(k(start=298.15*U.K))
          if internalHeat and diabatic constrainedby
          Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Heat",
            __Dymola_descriptionLabel=true,
            enable=(thermoOpt == 2 or thermoOpt == 3) and diabatic),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={50,20})));

        FCSys.Connectors.FaceGeneric negative(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Single-species connector for material, linear momentum, and heat"
          annotation (Placement(transformation(extent={{-110,-60},{-90,-40}}),
              iconTransformation(extent={{-110,-10},{-90,10}})));
        FCSys.Connectors.FaceGeneric positive(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Single-species connector for material, linear momentum, and heat"
          annotation (Placement(transformation(extent={{90,-60},{110,-40}}),
              iconTransformation(extent={{90,-10},{110,10}})));
        FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalHeat), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,70}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        Connectors.RealInputInternal u_N if not internalMaterial and open
          "External material signal" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-50,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-60,40})));
        Connectors.RealInputInternal u_Phi_1 if not internalLin1 and slip1
          "External signal for 1st transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Phi_2 if not internalLin2 and slip2
          "External signal for 2nd transverse linear momentum" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,40})));
        Connectors.RealInputInternal u_Q if not internalHeat and diabatic
          "External signal for heat" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={70,40}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,40})));

      equation
        // Material
        connect(materialBC.negative, negative.material) annotation (Line(
            points={{-70,-20},{-80,-20},{-80,-50},{-100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(materialBC.positive, positive.material) annotation (Line(
            points={{-50,-20},{80,-20},{80,-50},{100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(materialSpec.y, materialBC.u) annotation (Line(
            points={{-70,9},{-70,0},{-60,0},{-60,-16}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_N, materialBC.u) annotation (Line(
            points={{-50,40},{-50,0},{-60,0},{-60,-16}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.N, u_N) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-50,50},{-50,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 1st transverse linear momentum
        connect(lin1BC.negative, negative.mechanical1) annotation (Line(
            points={{-30,-30},{-80,-30},{-80,-50},{-100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin1BC.positive, positive.mechanical1) annotation (Line(
            points={{-10,-30},{80,-30},{80,-50},{100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin1Spec.y, lin1BC.u) annotation (Line(
            points={{-30,9},{-30,0},{-20,0},{-20,-26}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_1, lin1BC.u) annotation (Line(
            points={{-10,40},{-10,0},{-20,0},{-20,-26}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_1, u_Phi_1) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{-10,50},{-10,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // 2nd transverse linear momentum
        connect(lin2BC.negative, negative.mechanical2) annotation (Line(
            points={{10,-40},{-80,-40},{-80,-50},{-100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin2BC.positive, positive.mechanical2) annotation (Line(
            points={{30,-40},{80,-40},{80,-50},{100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(lin2Spec.y, lin2BC.u) annotation (Line(
            points={{10,9},{10,0},{20,0},{20,-36}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Phi_2, lin2BC.u) annotation (Line(
            points={{30,40},{30,0},{20,0},{20,-36}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Phi_2, u_Phi_2) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{30,50},{30,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Heat
        connect(heatBC.negative, negative.thermal) annotation (Line(
            points={{50,-50},{-100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(heatBC.positive, positive.thermal) annotation (Line(
            points={{70,-50},{100,-50}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(heatSpec.y, heatBC.u) annotation (Line(
            points={{50,9},{50,0},{60,0},{60,-46}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u_Q, heatBC.u) annotation (Line(
            points={{70,40},{70,0},{60,0},{60,-46}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(u.Q, u_Q) annotation (Line(
            points={{5.55112e-16,70},{5.55112e-16,50},{70,50},{70,40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        annotation (
          defaultComponentName="speciesFaceBC",
          Diagram(graphics),
          Icon(graphics));
      end Species;

      package Material "BCs for material"
        extends Modelica.Icons.Package;

        model PotentialElectrochemical "Prescribed electrochemical potential"

          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.PotentialElectrochemicalDifference,
              redeclare Connectors.RealInput u(final unit="l2.m/(N.T2)"));

        equation
          negative.mu - positive.mu = u "Condition";
          0 = negative.Ndot + positive.Ndot
            "Material rate balance (no storage)";
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="materialBC",
            Documentation(info="<html><p>Pressure can be
  calculated from other properties (e.g., temperature, specific volume, specific enthalpy, specific heat, or Gibbs potential)
  using functions in the
  <a href=\"modelica://FCSys.Connectors.Characteristic\">Characteristics</a> package.
  </p></html>"));
        end PotentialElectrochemical;

        model Current "Prescribed current"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Current,
              redeclare Connectors.RealInput u(final unit="N/T"));

        equation
          negative.Ndot = u "Condition";
          0 = negative.Ndot + positive.Ndot
            "Material rate balance (no storage)";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="materialBC");
        end Current;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a material BC"
            extends FCSys.BaseClasses.Icons.BCs.Double;
            constant BCType bCType "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.MaterialTransport negative
              "Material connector for the negative face" annotation (Placement(
                  transformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.MaterialTransport positive
              "Material connector for the positive face" annotation (Placement(
                  transformation(extent={{90,-10},{110,10}})));
            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="materialBC",
              Diagram(graphics));
          end PartialBC;

          type BCType = enumeration(
              PotentialElectrochemicalDifference
                "Electrochemical potential difference",
              Current "Current") "Types of BCs";
        end BaseClasses;
      end Material;

      package Momentum "BCs for linear momentum"
        extends Modelica.Icons.Package;
        model Velocity "Prescribed relative velocity"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.VelocityDifference,
              redeclare Connectors.RealInput u(final unit="l/T"));
        equation
          negative.phi - positive.phi = u "Condition";
          0 = negative.mPhidot + positive.mPhidot
            "Linear momentum rate balance (no storage)";

          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end Velocity;

        model Force "Prescribed total force with uniform velocity"
          extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force,
              redeclare Connectors.RealInput u(final unit="l.m/T2"));

        equation
          negative.mPhidot + positive.mPhidot = u "Condition";
          negative.phi = positive.phi "Uniform velocity";

          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="momBC");
        end Force;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a BC for linear momentum"
            extends FCSys.BaseClasses.Icons.BCs.Double;
            constant
              FCSys.BCs.FaceDifferential.Species.Momentum.BaseClasses.BCType
              bCType "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.MechanicalTransport negative
              "Linear momentum connector for the negative face" annotation (
                Placement(transformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.MechanicalTransport positive
              "Linear momentum connector for the positive face" annotation (
                Placement(transformation(extent={{90,-10},{110,10}})));
            annotation (defaultComponentPrefixes="replaceable",
                defaultComponentName="momBC");
          end PartialBC;

          type BCType = enumeration(
              VelocityDifference "Relative velocity",
              Force "Force") "Types of BCs";
        end BaseClasses;
      end Momentum;

      package Heat "BCs for heat"
        extends Modelica.Icons.Package;

        model Temperature "Prescribed temperature difference"
          extends FCSys.BCs.FaceDifferential.Species.Heat.BaseClasses.PartialBC(
              final bCType=BaseClasses.BCType.TemperatureDifference, redeclare
              Connectors.RealInput u(final unit="l2.m/(N.T2)", displayUnit="K"));
        equation
          negative.T - positive.T = u "Condition";
          0 = negative.Qdot + positive.Qdot "Energy rate balance (no storage)";
          annotation (defaultComponentPrefixes="replaceable",
              defaultComponentName="heatBC");
        end Temperature;

        model HeatFlowRate "Prescribed heat flow rate"
          extends FCSys.BCs.FaceDifferential.Species.Heat.BaseClasses.PartialBC(
              final bCType=BaseClasses.BCType.HeatFlowRate, redeclare
              Connectors.RealInput u(final unit="l2.m/T3"));

        equation
          negative.Qdot = u "Condition";
          0 = negative.Qdot + positive.Qdot "Energy rate balance (no storage)";

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="heatBC",
            Diagram(graphics));
        end HeatFlowRate;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;
          partial model PartialBC "Partial model for a BC for heat"
            extends FCSys.BaseClasses.Icons.BCs.Double;
            constant BCType bCType "Type of BC";
            // Note:  This is included so that the type of BC is recorded with the
            // results.
            Connectors.RealInput u "Value of BC" annotation (Placement(
                  transformation(
                  extent={{-10,-10},{10,10}},
                  rotation=270,
                  origin={0,40})));

            FCSys.Connectors.Thermal negative
              "Heat connector for the negative face" annotation (Placement(
                  transformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.Thermal positive
              "Heat connector for the positive face" annotation (Placement(
                  transformation(extent={{90,-10},{110,10}})));
          equation

            annotation (
              defaultComponentPrefixes="replaceable",
              defaultComponentName="heatBC",
              Diagram(graphics));
          end PartialBC;

          type BCType = enumeration(
              TemperatureDifference "Temperature difference",
              HeatFlowRate "Heat flow rate") "Types of BCs";
        end BaseClasses;
      end Heat;
    end Species;
    annotation (Documentation(info="<html><p>The hierarchy of these
 boundary condition models is similar to that of the models in the
 <a href=\"modelica://FCSys.BCs.Face\">Face boundary conditions</a> package.
 For more information, please see the documentation in that package.</p></html>"));
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
      Icon(graphics={
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
            smooth=Smooth.Bezier)}),
      Diagram(graphics));
  end Router;

  record Defaults "Global defaults for a model"
    extends FCSys.BaseClasses.Icons.Names.Top3;

    // Store the values of the base constants and units.
    final constant U.Bases.Base base=U.base "Base constants and units";

    parameter Boolean analysis=true "Include optional variables for analysis"
      annotation (choices(__Dymola_checkBox=true));

    parameter Q.PressureAbsolute p(nominal=1*U.atm) = 1*U.atm "Pressure";
    parameter Q.TemperatureAbsolute T(nominal=298.15*U.K) = 298.15*U.K
      "Temperature";
    parameter Q.NumberAbsolute RH(displayUnit="%") = 1 "Relative humidity";
    parameter Q.NumberAbsolute x_O2_dry(
      final max=1,
      displayUnit="%") = 0.208
      "<html>Dry gas O<sub>2</sub> fraction (<i>y</i><sub>O2 dry</sub>)</html>";
    // Value from http://en.wikipedia.org/wiki/Oxygen

    final parameter Q.NumberAbsolute x_H2O(
      final max=1,
      displayUnit="%") = 0.2
      "<html>Gas H<sub>2</sub>O fraction (<i>y</i><sub>H2O</sub>)</html>";
    // TODO:  Cast this in terms of relative humidity.

    annotation (
      defaultComponentPrefixes="inner",
      missingInnerMessage="Your model is using an outer \"defaults\" record, but an inner \"defaults\" record is not defined.
For simulation, specify global default settings by dragging FCSys.BCs.Defaults into your model.
The default global default settings will be used for the current simulation.",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            extent={{-21.2132,21.2132},{21.2131,-21.2131}},
            lineColor={127,127,127},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={225,225,225},
            origin={30,30},
            rotation=135),
          Polygon(
            points={{-60,60},{-60,-100},{60,-100},{60,30},{30,30},{30,60},{-60,
                60}},
            lineColor={127,127,127},
            smooth=Smooth.None,
            fillColor={245,245,245},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-60,60},{-60,-100},{60,-100},{60,30},{30,60},{-60,60}},
            lineColor={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-40,4},{40,4}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-40,-20},{40,-20}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-40,-44},{40,-44}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-40,-68},{40,-68}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-40,28},{26,28}},
            color={127,127,127},
            smooth=Smooth.None)}),
      Diagram(graphics));
  end Defaults;

  model ClosedVolume
    "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
    // Copied from Modelica.Fluid.

    import Modelica.Constants.pi;
    // Mass and energy balance, ports
    extends BaseClasses.PartialLumpedVessel(
      final fluidVolume=V,
      vesselArea=pi*(3/4*V)^(2/3),
      heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));
    parameter SI.Volume V "Volume";
  equation
    Wb_flow = 0;
    for i in 1:nPorts loop
      vessel_ps_static[i] = medium.p;
    end for;
    annotation (
      defaultComponentName="volume",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={170,213,255}), Text(
            extent={{-150,12},{150,-18}},
            lineColor={0,0,0},
            textString="V=%V")}),
      Documentation(info="<html>
<p>
Ideally mixed volume of constant size with two fluid ports and one medium model.
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>.
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected.
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</p>
<p>
If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between volume and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>.
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics));
  end ClosedVolume;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    partial model PartialLumpedVessel
      "Lumped volume with a vector of fluid ports and replaceable heat transfer model"
      // Copied and modified from Modelica.Fluid.
      import FCSys.Units;

      extends Modelica.Fluid.Interfaces.PartialLumpedVolume;
      // Port definitions
      parameter Integer nPorts=0 "Number of ports" annotation (Evaluate=true,
          Dialog(
          connectorSizing=true,
          tab="General",
          group="Ports"));
      Modelica.Fluid.Vessels.BaseClasses.VesselFluidPorts_b ports[nPorts](
          redeclare each package Medium = Medium) "Fluid inlets and outlets"
        annotation (Placement(transformation(extent={{-40,-10},{40,10}}, origin
              ={0,-100})));
      // Port properties
      parameter Boolean use_portsData=true
        "= false to neglect pressure loss and kinetic energy"
        annotation (Evaluate=true,Dialog(tab="General",group="Ports"));
      parameter Modelica.Fluid.Vessels.BaseClasses.VesselPortsData[nPorts]
        portsData if use_portsData "Data of inlet/outlet ports" annotation (
          Dialog(
          tab="General",
          group="Ports",
          enable=use_portsData));
      parameter SI.MassFlowRate m_flow_small(min=0) = system.m_flow_small
        "Regularization range at zero mass flow rate" annotation (Dialog(
          tab="Advanced",
          group="Port properties",
          enable=stiffCharacteristicForEmptyPort));
      /*
  parameter Medium.AbsolutePressure dp_small = system.dp_small
    "Turbulent flow if |dp| >= dp_small (regularization of zero flow)"
    annotation(Dialog(tab="Advanced",group="Ports"));
*/
      Medium.EnthalpyFlowRate ports_H_flow[nPorts];
      Medium.MassFlowRate ports_mXi_flow[nPorts, Medium.nXi];
      Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
        "Substance mass flows through ports";
      Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts, Medium.nC];
      Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
        "Trace substance mass flows through ports";
      // Heat transfer through boundary
      parameter Boolean use_HeatTransfer=false
        "= true to use the HeatTransfer model"
        annotation (Dialog(tab="Assumptions", group="Heat transfer"));
      replaceable model HeatTransfer =
          Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
        constrainedby
        Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
        "Wall heat transfer" annotation (Dialog(
          tab="Assumptions",
          group="Heat transfer",
          enable=use_HeatTransfer), choicesAllMatching=true);
      HeatTransfer heatTransfer(
        redeclare final package Medium = Medium,
        final n=1,
        final states={medium.state},
        final use_k=use_HeatTransfer) annotation (Placement(transformation(
            extent={{-10,-10},{30,30}},
            rotation=90,
            origin={-50,-10})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if
        use_HeatTransfer
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      // Conservation of kinetic energy
      Medium.Density[nPorts] portDensities
        "densities of the fluid at the device boundary";
      SI.Velocity[nPorts] portVelocities
        "velocities of fluid flow at device boundary";
      SI.EnergyFlowRate[nPorts] ports_E_flow
        "flow of kinetic and potential energy at device boundary";
      // Note:  should use fluidLevel_start - portsData.height
      Real[nPorts] s(each start=fluidLevel_max)
        "curve parameters for port flows vs. port pressures; for further details see, Modelica Tutorial: Ideal switching devices";
      Real[nPorts] ports_penetration
        "penetration of port with fluid, depending on fluid level and port diameter";
      // treatment of pressure losses at ports
      SI.Area[nPorts] portAreas={Modelica.Constants.pi/4*portsData_diameter[i]^
          2 for i in 1:nPorts};
      Medium.AbsolutePressure[nPorts] vessel_ps_static
        "static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

      Connectors.FaceGeneric face(
        final thermoOpt=ThermoOpt.OpenDiabatic,
        final slip1=false,
        final slip2=false)
        "Connection to a face of a FCSys.Subregions.Species model"
        annotation (Placement(transformation(extent={{90,-10},{110,10}})));

    protected
      input SI.Height fluidLevel=0
        "level of fluid in the vessel for treating heights of ports";
      parameter SI.Height fluidLevel_max=1
        "maximum level of fluid in the vessel";
      parameter SI.Area vesselArea=Modelica.Constants.inf
        "Area of the vessel used to relate to cross flow area of ports";
      // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
      // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
      // providing portsData_diameter and portsData_height, independent of the use_portsData setting.
      // Note:  this moreover serves as work-around if a tool doesn't support a zero sized portsData record.
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter_internal=
          portsData.diameter if use_portsData and nPorts > 0;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal=
          portsData.height if use_portsData and nPorts > 0;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal=
          portsData.zeta_in if use_portsData and nPorts > 0;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out_internal=
          portsData.zeta_out if use_portsData and nPorts > 0;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
      Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;
    equation
      // Added for FCSys:
      face.material.mu = Medium.specificGibbsEnergy(medium.state)*medium.MM*
        Units.J/Units.mol;
      face.thermal.T = medium.T*Units.K;

      mb_flow = sum(ports.m_flow) + medium.MM*face.material.Ndot/Units.kat
        "Changed for FCSys";
      mbXi_flow = sum_ports_mXi_flow;
      mbC_flow = sum_ports_mC_flow;
      Hb_flow = sum(ports_H_flow) + sum(ports_E_flow) + medium.h*medium.MM*face.material.Ndot
        /Units.kat "Changed for FCSys";
      Qb_flow = heatTransfer.Q_flows[1];
      // Only one connection allowed to a port to avoid unwanted ideal mixing
      for i in 1:nPorts loop
        assert(cardinality(ports[i]) <= 1, "
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeler. Increase nPorts to add an additional port.
");
      end for;
      // Check for correct solution
      assert(fluidLevel <= fluidLevel_max,
        "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(
        fluidLevel) + ")");
      assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(
        fluidLevel) + ") is below zero meaning that the solution failed.");
      // Boundary conditions
      // treatment of conditional portsData
      connect(portsData_diameter, portsData_diameter_internal);
      connect(portsData_height, portsData_height_internal);
      connect(portsData_zeta_in, portsData_zeta_in_internal);
      connect(portsData_zeta_out, portsData_zeta_out_internal);
      if not use_portsData then
        portsData_diameter = zeros(nPorts);
        portsData_height = zeros(nPorts);
        portsData_zeta_in = zeros(nPorts);
        portsData_zeta_out = zeros(nPorts);
      end if;
      // actual definition of port variables
      for i in 1:nPorts loop
        if use_portsData then
          // dp = 0.5*zeta*d*v*|v|
          // Note:  assume vessel_ps_static for portDensities to avoid algebraic loops for ports.p
          portDensities[i] = noEvent(Medium.density(Medium.setState_phX(
                vessel_ps_static[i],
                actualStream(ports[i].h_outflow),
                actualStream(ports[i].Xi_outflow))));
          portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/
            portDensities[i]);
          // Note:  the penetration should not go too close to zero as this would prevent a vessel from running empty
          ports_penetration[i] = Modelica.Fluid.Utilities.regStep(
                fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i],
                1,
                1e-3,
                0.1*portsData_diameter[i]);
        else
          // an infinite port diameter is assumed
          portDensities[i] = medium.d;
          portVelocities[i] = 0;
          ports_penetration[i] = 1;
        end if;
        // fluid flow through ports
        if fluidLevel >= portsData_height[i] then
          // regular operation: fluidLevel is above ports[i]
          // Note:  >= covers default values of zero as well
          if use_portsData then
            /* Without regularization
        ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2
                      * noEvent(if ports[i].m_flow>0 then zeta_in[i]/portDensities[i] else -zeta_out[i]/medium.d);
        */
            ports[i].p = vessel_ps_static[i] + (0.5/portAreas[i]^2*
              Modelica.Fluid.Utilities.regSquare2(
                  ports[i].m_flow,
                  m_flow_small,
                  (portsData_zeta_in[i] - 1 + portAreas[i]^2/vesselArea^2)/
                portDensities[i]*ports_penetration[i],
                  (portsData_zeta_out[i] + 1 - portAreas[i]^2/vesselArea^2)/
                medium.d/ports_penetration[i]));
            /*
        // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
        ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                     2*portDensities[i]/portsData_zeta_in[i],
                                     2*medium.d/portsData_zeta_out[i]));
        */
          else
            ports[i].p = vessel_ps_static[i];
          end if;
          s[i] = fluidLevel - portsData_height[i];
        elseif s[i] > 0 or portsData_height[i] >= fluidLevel_max then
          // ports[i] is above fluidLevel and has inflow
          ports[i].p = vessel_ps_static[i];
          s[i] = ports[i].m_flow;
        else
          // ports[i] is above fluidLevel, preventing outflow
          ports[i].m_flow = 0;
          s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(
            portsData_height[i] - fluidLevel);
        end if;
        ports[i].h_outflow = medium.h;
        ports[i].Xi_outflow = medium.Xi;
        ports[i].C_outflow = C;
        ports_H_flow[i] = ports[i].m_flow*actualStream(ports[i].h_outflow)
          "Enthalpy flow";
        ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[
          i] + system.g*portsData_height[i])
          "Flow of kinetic and potential energy";
        ports_mXi_flow[i, :] = ports[i].m_flow*actualStream(ports[i].Xi_outflow)
          "Component mass flow";
        ports_mC_flow[i, :] = ports[i].m_flow*actualStream(ports[i].C_outflow)
          "Trace substance mass flow";
      end for;
      for i in 1:Medium.nXi loop
        sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:, i]);
      end for;
      for i in 1:Medium.nC loop
        sum_ports_mC_flow[i] = sum(ports_mC_flow[:, i]);
      end for;
      connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
          points={{-100,5.55112e-16},{-87,5.55112e-16},{-87,2.22045e-15},{-74,
              2.22045e-15}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (
        Documentation(info="<html>
<p>
This base class extends PartialLumpedVolume with a vector of fluid ports and a replaceable wall HeatTransfer model.
<p>
The following modeling assumption are made:
<ul>
<li>homogeneous medium, i.e., phase separation is not taken into account,</li>
<li>no kinetic energy in the fluid, i.e., kinetic energy dissipates into the internal energy,</li>
<li>pressure loss definitions at vessel ports assume incompressible fluid,</li>
<li>outflow of ambient media is prevented at each port assuming check valve behavior.
    If <code> fluidlevel &lt; portsData_height[i] </code>and &nbsp; <code> ports[i].p &lt; vessel_ps_static[i]</code> massflow at the port is set to 0.</li>
</ul>
</p>
Each port has a (hydraulic) diameter and a height above the bottom of the vessel, which can be configured using the &nbsp;<b><code>portsData</code></b> record.
Alternatively the impact of port geometries can be neglected with <code>use_portsData=false</code>. This might be useful for early
design studies. Note that this means to assume an infinite port diameter at the bottom of the vessel.
Pressure drops and heights of the ports as well as kinetic and potential energy fluid entering or leaving the vessel are neglected then.
<p>
The following variables need to be defined by an extending model:
<ul>
<li><code>input fluidVolume</code>, the volume of the fluid in the vessel,</li>
<li><code>vessel_ps_static[nPorts]</code>, the static pressures inside the vessel at the height of the corresponding ports, at zero flow velocity, and</li>
<li><code>Wb_flow</code>, work term of the energy balance, e.g., p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
An extending model should define:
<ul>
<li><code>parameter vesselArea</code> (default: Modelica.Constants.inf m2), the area of the vessel, to be related to cross flow areas of the ports for the consideration of dynamic pressure effects.</li>
</ul>
Optionally the fluid level may vary in the vessel, which effects the flow through the ports at configurable <code>portsData_height[nPorts]</code>.
This is why an extending model with varying fluid level needs to define:
<ul>
<li><code>input fluidLevel (default: 0m)</code>, the level the fluid in the vessel, and</li>
<li><code>parameter fluidLevel_max (default: 1m)</code>, the maximum level that must not be exceeded. Ports at or above fluidLevel_max can only receive inflow.</li>
</ul>
An extending model should not access the <code>portsData</code> record defined in the configuration dialog,
as an access to <code>portsData</code> may fail for <code>use_portsData=false</code> or <code>nPorts=0</code>.
Instead the predefined variables
<ul>
<li><code>portsData_diameter[nPorts]</code></li>,
<li><code>portsData_height[nPorts]</code></li>,
<li><code>portsData_zeta_in[nPorts]</code></li>, and
<li><code>portsData_zeta_out[nPorts]</code></li>
</ul>
should be used if these values are needed.
</p>
</html>", revisions="<html>
<ul>
<li><i>Jan. 2009</i> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic and potential energy of fluid entering or leaving in energy balance</li>
   </ul>
</li>
<li><i>Dec. 2008</i> by R&uuml;diger Franke: derived from OpenTank, in order to make general use of configurable port diameters</i>
</ul>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
                100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Text(
              extent={{-150,110},{150,150}},
              textString="%name",
              lineColor={0,0,255})}));
    end PartialLumpedVessel;

    block RealFunction
      "<html>Set output signal according to a <code>Real</code> function of an input</html>"

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
