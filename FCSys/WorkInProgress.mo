within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;
  extends FCSys.Icons.PackageUnderConstruction;

  model ConditionsAdaptersPhasesIonomer
    "<html>Adapter for ionomer between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
    extends Conditions.Adapters.Phases.BaseClasses.PartialPhase;
    extends Modelica.Icons.UnderConstruction;
    Conditions.Adapters.Species.Solid 'SO3-';
    FCSys.WorkInProgress.ConditionsAdaptersSpeciesFluid 'H+'(redeclare package
        Data = FCSys.Characteristics.'H+'.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2)
      annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
    // **Use model for H instead.
    Conditions.Adapters.Species.FluidNeutral H2O(redeclare package Data =
          Characteristics.H2O.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2O)
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Conditions.Adapters.Junctions.Junction2 junction2
      annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
        Medium = Medium) "Modelica fluid port" annotation (Placement(
          transformation(extent={{70,-50},{90,-30}}), iconTransformation(extent
            ={{70,-50},{90,-30}})));
    Modelica.Electrical.Analog.Interfaces.NegativePin pin
      "Modelica electrical pin" annotation (Placement(transformation(extent={{
              70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

  equation
    connect('SO3-'.face.thermal, face.'SO3-'.thermal) annotation (Line(
        points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.face.normal, face.'H+'.normal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.pin, pin) annotation (Line(
        points={{8,-16},{60,-16},{60,40},{80,40}},
        color={0,0,255},
        smooth=Smooth.None));
    connect('H+'.heatPort, heatPort) annotation (Line(
        points={{8,-20},{40,-20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect('SO3-'.heatPort, heatPort) annotation (Line(
        points={{8,20},{40,20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(junction2.mixturePort, fluidPort) annotation (Line(
        points={{58,-40},{80,-40}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (Placement(transformation(extent={{-10,10},{10,30}})), Icon(
          graphics={Line(
              points={{0,60},{0,-60}},
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
              smooth=Smooth.None),Line(
              points={{0,-40},{80,-40}},
              color={0,127,255},
              smooth=Smooth.None)}));
  end ConditionsAdaptersPhasesIonomer;

  model CellModelica
    "<html>Cell interfaced to components from the <a href=\"modelica://Modelica\">Modelica</a> package</html>"
    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    Assemblies.Cells.Examples.Cell cell(anFP(redeclare
          FCSys.Subregions.Subregion subregions(
          each final inclX=true,
          each inclY=true,
          each graphite('incle-'=true, 'e-'(perfectMaterialDiff={{{{true,false}}}})),

          each gas(inclH2=true, inclH2O=true))))
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    inner FCSys.Conditions.Environment environment(analysis=false)
      annotation (Placement(transformation(extent={{40,60},{60,80}})));
    Conditions.Adapters.Phases.Graphite caModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1]) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    Conditions.Adapters.Phases.Graphite anModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1])
      annotation (Placement(transformation(extent={{-20,-10},{-40,10}})));
    FCSys.WorkInProgress.TanConduct tanConduct
      annotation (Placement(transformation(extent={{10,40},{-10,60}})));
    Modelica.Blocks.Sources.Ramp loadSweep(duration=1000)
      "This is the arctangent of conductance."
      annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

  equation
    connect(anModelicaAdapt.normal, cell.anFPX[1, 1]) annotation (Line(
        points={{-20,6.10623e-16},{-16,6.10623e-16},{-16,5.55112e-16},{-10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(caModelicaAdapt.normal, cell.caFPX[1, 1]) annotation (Line(
        points={{20,6.10623e-16},{16,6.10623e-16},{16,5.55112e-16},{10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(loadSweep.y, tanConduct.atanGstar) annotation (Line(
        points={{-19,70},{-6.66134e-16,70},{-6.66134e-16,61}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tanConduct.p, caModelicaAdapt.pin) annotation (Line(
        points={{10,50},{60,50},{60,4},{40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tanConduct.n, anModelicaAdapt.pin) annotation (Line(
        points={{-10,50},{-60,50},{-60,4},{-40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (experiment(StopTime=1000));
  end CellModelica;

  function plot "Create plots using FCRes"
    extends Modelica.Icons.Function;
    extends Modelica.Icons.UnderConstruction;

  algorithm
    Modelica.Utilities.System.command("loadres");

  end plot;

  model EISPlaceholder
    "Placeholder model for electrochemical-impedance spectroscopy"
    extends FCSys.Icons.Blocks.Continuous;

    parameter Modelica.SIunits.CurrentDensity zJ_large_SI=0.01
      "Large-signal current density in SI base units";
    Modelica.Electrical.Analog.Basic.Resistor resistor2(R=0.1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,-10})));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=1e-3) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,-10})));
    Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=1)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-80,-10})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={70,-10})));
    Modelica.Electrical.Analog.Sources.ConstantCurrent constantCurrent(I=
          zJ_large_SI)
      annotation (Placement(transformation(extent={{-60,20},{-40,0}})));
    Modelica.Electrical.Analog.Sources.SignalCurrent signalCurrent
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
    Connectors.RealInput zJ_small_SI
      "Small-signal cell current density in SI base units" annotation (
        Placement(transformation(extent={{-110,40},{-90,60}}),
          iconTransformation(extent={{-120,-10},{-100,10}})));
    Connectors.RealOutput w_V "Cell potential in volts" annotation (Placement(
          transformation(extent={{90,-20},{110,0}}), iconTransformation(extent=
              {{100,-10},{120,10}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor1(R=0.1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,10})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{60,-60},{80,-40}})));
  equation
    connect(zJ_small_SI, signalCurrent.i) annotation (Line(
        points={{-100,50},{-50,50},{-50,37}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(signalCurrent.n, constantCurrent.n) annotation (Line(
        points={{-40,30},{-30,30},{-30,10},{-40,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(signalCurrent.p, constantCurrent.p) annotation (Line(
        points={{-60,30},{-70,30},{-70,10},{-60,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantVoltage.p, constantCurrent.p) annotation (Line(
        points={{-80,5.55112e-16},{-80,10},{-60,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(voltageSensor.v, w_V) annotation (Line(
        points={{80,-10},{90,-10},{90,-10},{100,-10}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(resistor2.p, capacitor.p) annotation (Line(
        points={{10,5.55112e-16},{10,10},{40,10},{40,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantVoltage.n, resistor2.n) annotation (Line(
        points={{-80,-20},{-80,-30},{10,-30},{10,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.n, resistor2.n) annotation (Line(
        points={{40,-20},{40,-30},{10,-30},{10,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.n, voltageSensor.n) annotation (Line(
        points={{40,-20},{40,-30},{70,-30},{70,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.p, voltageSensor.p) annotation (Line(
        points={{40,5.55112e-16},{40,10},{70,10},{70,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantCurrent.n, resistor1.p) annotation (Line(
        points={{-40,10},{-20,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistor1.n, resistor2.p) annotation (Line(
        points={{5.55112e-16,10},{10,10},{10,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(ground.p, voltageSensor.n) annotation (Line(
        points={{70,-40},{70,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (Diagram(graphics), Icon(graphics));
  end EISPlaceholder;

  model ChargeLayerStoich
    "<html>Adapter between the <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> and <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors</html>"

    import FCSys.Utilities.inSign;
    extends FCSys.Icons.Names.Top2;

    // Geometry
    parameter Q.Area A "Cross-sectional area of the face" annotation (Dialog(
          group="Geometry", __Dymola_label="<html><i>A</i></html>"));
    parameter Axis axis "Axis of the electrochemical reaction"
      annotation (Dialog(group="Geometry"));
    parameter Side side "Side of the face w.r.t., the reaction"
      annotation (Dialog(group="Geometry"));
    parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (Dialog(group="Geometry"));
    parameter Integer n "Stoichiometric coefficient"
      annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));

    replaceable package Data = Characteristics.BaseClasses.Characteristic
      constrainedby Characteristics.BaseClasses.Characteristic
      "Characteristic data" annotation (
      Dialog(group="Material properties"),
      choicesAllMatching=true,
      __Dymola_choicesFromPackage=true);

    FCSys.Connectors.Reaction reaction(final n_trans=n_trans)
      "Connector for an electrochemical reaction" annotation (
        __Dymola_choicesAllMatching=true, Placement(transformation(extent={{30,
              -10},{50,10}}), iconTransformation(extent={{30,-10},{50,10}})));

    FCSys.Conditions.Adapters.ChargeLayer chargeLayer(
      final axis=axis,
      final side=side,
      final cartTrans=cartTrans,
      redeclare final package Data = Data,
      final A=A)
      annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
    FCSys.Conditions.Adapters.Stoich stoich(
      final n_trans=n_trans,
      final n=n,
      final m=Data.m)
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));

  protected
    final parameter Integer n_trans=size(cartTrans, 1)
      "Number of components of translational momentum";

  public
    Connectors.Face substrate "Boundary of the substrate" annotation (Placement(
          transformation(extent={{-50,-30},{-30,-10}}), iconTransformation(
            extent={{-50,-30},{-30,-10}})));
  public
    Connectors.Face carrier "Boundary of the charge carrier" annotation (
        Placement(transformation(extent={{-50,10},{-30,30}}),
          iconTransformation(extent={{-50,10},{-30,30}})));
  equation

    connect(stoich.reaction, reaction) annotation (Line(
        points={{14,6.10623e-16},{40,5.55112e-16}},
        color={255,195,38},
        smooth=Smooth.None));

    connect(carrier, chargeLayer.carrier) annotation (Line(
        points={{-40,20},{-30,20},{-30,2},{-14,2}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(chargeLayer.substrate, substrate) annotation (Line(
        points={{-14,-2},{-30,-2},{-30,-20},{-40,-20}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(chargeLayer.electrical, stoich.chemical) annotation (Line(
        points={{-6,6.10623e-16},{6,6.10623e-16},{6,6.10623e-16}},
        color={255,195,38},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html><p>This model provides the stoichiometry-weighted electrical contribution of a charge
    carrier to an electrochemical reaction.  The <code>carrier</code> and <code>substrate</code>
    connectors represent the interface with the charge carrier species (e.g., e<sup>-</sup> or H<sup>+</sup>)
    and the solid species (e.g., C<sup>+</sup>) at concentration at the far edge of the far edge of the 
    depletion region.  The concentration of the charge carrier there results in the Butler-Volmer equation 
    in conjunction with the material transport of the charge carrier species.</p>  
    
    <p>Assumptions:<ol>
    <li>There is no diffusion of material, transverse translational momentum, or energy across the face.</li>
    <li>The diffusive or non-equilibrium normal force is applied to the electrical part of the electrochemical
    potential.</li> 
    </ol></p>
    
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics),
      Icon(graphics={Line(
              points={{0,0},{30,0}},
              color={255,195,38},
              smooth=Smooth.None),Line(
              points={{-30,20},{0,20}},
              color={127,127,127},
              smooth=Smooth.None),Line(
              points={{0,-30},{0,30}},
              color={127,127,127},
              smooth=Smooth.None,
              thickness=0.5),Text(
              extent={{-100,-40},{100,-60}},
              lineColor={127,127,127},
              textString="%n"),Line(
              points={{-30,-20},{0,-20}},
              color={127,127,127},
              smooth=Smooth.None)}));
  end ChargeLayerStoich;

  model StoichMulti
    "<html>Adapter between <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connectors and a <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector</html>"

    extends FCSys.Icons.Names.Top1;

    parameter Integer n_trans(
      min=1,
      max=3) = 1 "Number of components of translational momentum" annotation (
        Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
    parameter Integer n_spec "Number of species"
      annotation (Dialog(connectorSizing=true));

    parameter Integer n[n_spec] "Stoichiometric coefficients"
      annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
    parameter Q.MassSpecific m[n_spec] "Specific masses" annotation (Dialog(
          group="Material properties", __Dymola_label="<html><i>m</i></html>"));

    Connectors.Chemical chemical[n_spec](each final n_trans=n_trans) if n_spec
       > 0 "Connector for species in a chemical reaction" annotation (Placement(
          transformation(extent={{-30,-10},{-10,10}}), iconTransformation(
            extent={{-50,-10},{-30,10}})));
    Connectors.Reaction reaction(final n_trans=n_trans) if n_spec > 0
      "Connector for an electrochemical reaction" annotation (Placement(
          transformation(extent={{10,-10},{30,10}}), iconTransformation(extent=
              {{30,-10},{50,10}})));

  protected
    Conditions.Adapters.Stoich stoich[n_spec](
      each final n_trans=n_trans,
      final n=n,
      final m=m) if n_spec > 0
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  equation
    connect(stoich.chemical, chemical) annotation (Line(
        points={{-4,6.10623e-16},{-12,6.10623e-16},{-12,5.55112e-16},{-20,
            5.55112e-16}},
        color={255,195,38},
        smooth=Smooth.None));

    for i in 1:n_spec loop
      connect(stoich[i].reaction, reaction) annotation (Line(
          points={{4,6.10623e-16},{10,6.10623e-16},{10,5.55112e-16},{20,
              5.55112e-16}},
          color={255,195,38},
          smooth=Smooth.None));

    end for;

    annotation (
      defaultComponentName="stoich",
      Documentation(info="<html><p>This model is used to add the stoichiometrically-weighted electrochemical potentials
    of species to the net electrochemical potential of a reaction.  The species are produced at the
    stoichiometrically-weighted rate of the reaction.</p>
    
    <p>See also <a href=\"modelica://FCSys.Conditions.Adapters.ChemicalReaction\">ChemicalReaction</a>.  
    For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Line(
              points={{-30,0},{30,0}},
              color={255,195,38},
              smooth=Smooth.None),Text(
              extent={{-100,-20},{100,-40}},
              lineColor={127,127,127},
              textString="%n"),Line(
              points={{0,-10},{0,10}},
              color={127,127,127},
              smooth=Smooth.None,
              thickness=0.5)}),
      Diagram(graphics));
  end StoichMulti;
  annotation (Commands(
      file="../../units.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate).",

      file="test/check.mos" "Check all of FCSys using Dymola's check function.",

      file="../../../LaTeX/Dissertation/Results/Cell/Simulation/sim.mos"));

end WorkInProgress;
