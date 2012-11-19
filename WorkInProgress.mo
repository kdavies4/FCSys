within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;
  extends FCSys.BaseClasses.Icons.PackageUnderConstruction;

  model Subregion "Subregion with all phases included"

    parameter Boolean inclReact=true "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(tab="Assumptions"),
      choices(__Dymola_checkBox=true));
    // Note:  This is listed above the extends clause so that it is
    // listed first in the parameter dialog.
    extends Subregions.BaseClasses.PartialSubregion;

    replaceable FCSys.Subregions.Phases.Gas gas(inclH2O=true, final inclLin={
          inclLinX,inclLinY,inclLinZ}) "Gas" annotation (Dialog(group="Phases"),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    // **Currently, both reactions must be included at this level.

    replaceable FCSys.Subregions.Phases.Graphite graphite(final inclLin={
          inclLinX,inclLinY,inclLinZ}) "Graphite" annotation (Dialog(group=
            "Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable FCSys.Subregions.Phases.Ionomer ionomer(final inclLin={inclLinX,
          inclLinY,inclLinZ},C19HF37O5S(initMethTemp=if graphite.inclC then
            InitMethScalar.None else InitMethScalar.Temperature,T(stateSelect=
              if graphite.inclC then StateSelect.never else StateSelect.prefer)))
      "Ionomer" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    // **temp stateselect, initmeth
    /*
  replaceable FCSys.Subregions.Phases.Liquid liquid(final inclLin={inclLinX,
        inclLinY,inclLinZ})  "Liquid" annotation (

    Dialog(group="Phases"),
    Placement(transformation(extent={{-10,-10},{10,10}})));
  */

    FCSys.Subregions.Reactions.Electrochemical HOR(final n_lin=n_lin, n_spec=3)
      if inclReact and (graphite.'incle-' and ionomer.'inclH+' and gas.inclH2
       and not (gas.inclO2 and gas.inclH2O)) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
    FCSys.Subregions.Reactions.Electrochemical ORR(final n_lin=n_lin, n_spec=4)
      if inclReact and (graphite.'incle-' and ionomer.'inclH+' and gas.inclO2
       and gas.inclH2O and not gas.inclH2) "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
    FCSys.WorkInProgress.Capacitor DL if graphite.'incle-' and ionomer.'inclH+'
      "Double layer capacitor" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-40,20})));

    /* **Note:  Multiple reactions cannot be included at once due to the following error:
Error: Failed to expand the variable HOR.chemical[1].mphi
Error: Failed to expand the variable HOR.chemical[2].mphi
Error: Failed to expand the variable ORR.chemical[1].mphi
Error: Failed to expand the variable ORR.chemical[2].mphi
  */

  protected
    Connectors.ChemicalBusInternal chemical
      "Internal connector to route electrochemical interactions"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}})));

  equation
    // Chemical interactions
    connect(HOR.chemical[1], chemical.'e-') annotation (Line(
        points={{-20,39.3333},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[2], chemical.'H+') annotation (Line(
        points={{-20,40},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[3], chemical.H2) annotation (Line(
        points={{-20,40.6667},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(ORR.chemical[1], chemical.'e-') annotation (Line(
        points={{-20,39.25},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[2], chemical.'H+') annotation (Line(
        points={{-20,39.75},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[3], chemical.O2) annotation (Line(
        points={{-20,40.25},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[4], chemical.H2O) annotation (Line(
        points={{-20,40.75},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(DL.ion1, chemical.'e-') annotation (Line(
        points={{-40,14},{-30,14},{-30,20},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}}));
    connect(DL.ion2, chemical.'H+') annotation (Line(
        points={{-40,26},{-30,26},{-30,20},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}}));

    // Gas
    connect(gas.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-24,-4.87687e-22},{-24,5.55112e-16},{-40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{5.55112e-16,
            40}},
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
    connect(graphite.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(graphite.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-24,-4.87687e-22},{-24,5.55112e-16},{-40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{5.55112e-16,
            40}},
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
    connect(ionomer.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(ionomer.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-24,-4.87687e-22},{-24,5.55112e-16},{-40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
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

    /*
  // Liquid
  connect(liquid.inert, volume.interaction) annotation (Line(
      points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
      color={0,180,0},
      smooth=Smooth.None,
      thickness=0.5));
  connect(liquid.xNegative, xNegative.liquid) annotation (Line(
      points={{-10,6.10623e-16},{-10,1.16573e-15},{-25,1.16573e-15},{-25,5.55112e-16},
          {-40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.xPositive, xPositive.liquid) annotation (Line(
      points={{10,6.10623e-16},{10,5.55112e-16},{40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.yNegative, yNegative.liquid) annotation (Line(
      points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.yPositive, yPositive.liquid) annotation (Line(
      points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{5.55112e-16,
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
      points={{-5,-5},{-20,-20}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  */

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html><p>Notes:<ul>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
  <li>H<sub>2</sub>O is included by default, since at least one species
  must be included.</li></ul>
</p><p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Icon(graphics),
      Diagram(graphics));
  end Subregion;

  model Capacitor "Model for an electrical capacitor"
    extends FCSys.BaseClasses.Icons.Names.Top2;

    parameter Integer n_lin=1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (Evaluate=true,HideResult=true);
    parameter Q.Amount DeltaQ_IC(min=-Modelica.Constants.inf) = 0
      "Initial charge difference" annotation (Dialog(group="Initialization"));
    input Q.Amount CT(nominal=300*U.micro*U.F*U.K) = 300*U.K*U.micro*U.F
      "Capacitance times temperature" annotation (Dialog);

    Integer z[2]=Chemistry.charge({ion1.formula,ion2.formula}) "Charge numbers";
    Q.Current Qdot[2] "Incoming electrical currents";
    Q.Amount DeltaQ(
      min=-Modelica.Constants.inf,
      stateSelect=StateSelect.prefer,
      final start=DeltaQ_IC,
      fixed=true) "Charge difference";
    //**stateSelect=StateSelect.prefer,

    Connectors.ChemicalInput ion1(each final n_lin=n_lin)
      "Connector for the 1st ion" annotation (Placement(transformation(extent={
              {-70,-10},{-50,10}}), iconTransformation(extent={{-70,-10},{-50,
              10}})));
    Connectors.ChemicalInput ion2(each final n_lin=n_lin)
      "Connector for the 2nd ion" annotation (Placement(transformation(extent={
              {50,-10},{70,10}}), iconTransformation(extent={{50,-10},{70,10}})));

  equation
    // Aliases
    Qdot = z .* {ion1.Ndot,ion2.Ndot};
    DeltaQ/2 = CT*(ion1.muPerT/z[1] - ion2.muPerT/z[2]);

    // Electrostatic energy storage
    der(DeltaQ)/U.s = Delta(Qdot);

    // No exchange of linear momentum
    ion1.mPhidot = ion2.mPhidot;

    // No thermal exchange
    ion1.Hdot = ion2.Hdot;

    // Conservation (no storage)
    0 = Sigma(Qdot) "Charge (net)";
    zeros(n_lin) = ion1.mPhidot + ion2.mPhidot "Linear momentum";
    0 = ion1.Hdot + ion2.Hdot "Energy (excluding electrostatic)";

    // This model is marked as structurally incomplete because
    // by default the chemical formulas are empty strings and
    // thus the charge numbers are zero.
    annotation (
      structurallyIncomplete=true,
      Documentation(info="<html>
    <p>**
    </p>
    </html>"),
      Icon(graphics={Line(
              points={{-60,0},{-10,0}},
              color={208,104,0},
              smooth=Smooth.None),Line(
              points={{10,0},{60,0}},
              color={208,104,0},
              smooth=Smooth.None),Line(
              points={{-10,36},{-10,-36}},
              color={208,104,0},
              smooth=Smooth.None),Line(
              points={{10,36},{10,-36}},
              color={208,104,0},
              smooth=Smooth.None)}),
      Diagram(graphics));
  end Capacitor;
  annotation (Commands(file="resources/scripts/units-values.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate)."));
end WorkInProgress;
