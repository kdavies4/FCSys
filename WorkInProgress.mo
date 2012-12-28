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

    replaceable FCSys.Subregions.Phases.Gas gas(inclH2O=true, final inclVel={
          inclVelX,inclVelY,inclVelZ}) "Gas" annotation (Dialog(group="Phases"),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    // **Currently, both reactions must be included at this level.

    replaceable FCSys.Subregions.Phases.Graphite graphite(final inclVel={
          inclVelX,inclVelY,inclVelZ}) "Graphite" annotation (Dialog(group=
            "Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable FCSys.Subregions.Phases.Ionomer ionomer(final inclVel={inclVelX,
          inclVelY,inclVelZ},C19HF37O5S(initMethTemp=if graphite.inclC then
            InitMethScalar.None else InitMethScalar.Temperature,T(stateSelect=
              if graphite.inclC then StateSelect.never else StateSelect.prefer)))
      "Ionomer" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    // **temp stateSelect, initMeth
    /*
  replaceable FCSys.Subregions.Phases.Liquid liquid(final inclVel={inclVelX,
        inclVelY,inclVelZ})  "Liquid" annotation (

    Dialog(group="Phases"),
    Placement(transformation(extent={{-10,-10},{10,10}})));
  */

    FCSys.Subregions.Reactions.Electrochemical HOR(final n_vel=n_vel, n_spec=3)
      if inclReact and (graphite.'incle-' and ionomer.'inclH+' and gas.inclH2
       and not (gas.inclO2 and gas.inclH2O)) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
    FCSys.Subregions.Reactions.Electrochemical ORR(final n_vel=n_vel, n_spec=4)
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

    parameter Integer n_vel=1
      "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
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

    Connectors.ChemicalInput ion1(each final n_vel=n_vel)
      "Connector for the 1st ion" annotation (Placement(transformation(extent={
              {-70,-10},{-50,10}}),iconTransformation(extent={{-70,-10},{-50,10}})));
    Connectors.ChemicalInput ion2(each final n_vel=n_vel)
      "Connector for the 2nd ion" annotation (Placement(transformation(extent={
              {50,-10},{70,10}}),iconTransformation(extent={{50,-10},{70,10}})));

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
    zeros(n_vel) = ion1.mPhidot + ion2.mPhidot "Linear momentum";
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

  model H2_O2_H2ODynamic
    "<html>Test the 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O reaction dynamically</html>"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    parameter Q.TemperatureAbsolute T=298.15*U.K "Temperature";
    parameter Q.Volume V=1*U.cm^3 "Volume";

    FCSys.Subregions.Reaction reaction(n_spec=3) "Chemical reaction"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    FCSys.WorkInProgress.SimpleSpecies species1(
      final V=V,
      final T=T,
      p_IC=0.4*U.atm,
      p(fixed=false),
      redeclare FCSys.Characteristics.H2.Gas Data(b_v=[1], specVolPow={-1,0}))
      "1st species" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-30,-24})));
    FCSys.WorkInProgress.SimpleSpecies species2(
      final V=V,
      final T=T,
      p_IC=1*U.atm,
      redeclare FCSys.Characteristics.O2.Gas Data(b_v=[1], specVolPow={-1,0}))
      "2nd species" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={0,-24})));
    FCSys.WorkInProgress.SimpleSpecies species3(
      final V=V,
      final T=T,
      p_IC=0.4*U.atm,
      redeclare FCSys.Characteristics.H2O.Gas Data(b_v=[1], specVolPow={-1,0}))
      "3rd species" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={30,-24})));

    inner FCSys.BCs.Defaults defaults(analysis=true)
      annotation (Placement(transformation(extent={{-90,70},{-70,90}})));

  equation
    connect(species1.chemical, reaction.chemical[1]) annotation (Line(
        points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,-0.666667}},

        color={170,0,0},
        smooth=Smooth.None));

    connect(species2.chemical, reaction.chemical[2]) annotation (Line(
        points={{9.89443e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
            5.55112e-16}},
        color={170,0,0},
        smooth=Smooth.None));

    connect(species3.chemical, reaction.chemical[3]) annotation (Line(
        points={{30,-20},{30,-10},{5.55112e-16,-10},{5.55112e-16,0.666667}},
        color={170,0,0},
        smooth=Smooth.None));

    annotation (
      Diagram(graphics),
      experiment(
        StopTime=0.01,
        NumberOfIntervals=5000,
        Tolerance=1e-06),
      experimentSetupOutput,
      Commands(file=
            "resources/scripts/Dymola/Subregions.Examples.H2_O2_H2ODynamic.mos"));
  end H2_O2_H2ODynamic;

  model H2O_H2ODynamic
    "<html>Test the H<sub>2</sub>O evaporation/condensation reaction dynamically</html>"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    parameter Q.TemperatureAbsolute T=298.15*U.K "Temperature";
    parameter Q.Volume V=1*U.cm^3 "Volume";

    FCSys.Subregions.Reaction reaction(n_spec=2) "Chemical reaction"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    FCSys.WorkInProgress.SimpleSpecies species1(
      final V=V,
      final T=T,
      redeclare FCSys.Characteristics.H2O.Liquid Data,
      p(fixed=false)) "1st species" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-30,-24})));
    FCSys.WorkInProgress.SimpleSpecies species2(
      final V=V,
      final T=T,
      p_IC=0.4*U.atm,
      redeclare FCSys.Characteristics.H2O.Gas Data(b_v=[1], specVolPow={-1,0}))
      "2nd species" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={0,-24})));

    inner FCSys.BCs.Defaults defaults(analysis=true)
      annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
  equation

    connect(species1.chemical, reaction.chemical[1]) annotation (Line(
        points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,-0.5}},
        color={170,0,0},
        smooth=Smooth.None));

    connect(species2.chemical, reaction.chemical[2]) annotation (Line(
        points={{9.89443e-16,-20},{0,-20},{0,0.5},{5.55112e-16,0.5}},
        color={170,0,0},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>Note that the pressure of liquid H<sub>2</sub>O shows as zero.
  Since it is assumed to be incompressible, its pressure cannot be determined by the specific volume and temperature.
  In that case, the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p\">p</a> function
  returns zero.</html>"),
      Diagram(graphics),
      experiment(
        StopTime=2,
        NumberOfIntervals=5000,
        Tolerance=1e-06),
      experimentSetupOutput,
      Commands(file=
            "resources/scripts/Dymola/Subregions.Examples.H2O_H2ODynamic.mos"));
  end H2O_H2ODynamic;

  model Cell "Test both half reactions of a cell"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    output Q.Number DeltamuPerT='e-_an'.chemical.muPerT - 'e-_ca'.chemical.muPerT
      "Voltage of cathode w.r.t. anode, divided by temperature";
    // The negative factor is due to the charge negative charge of electrons.

    FCSys.Subregions.Reaction HOR(n_spec=3) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

    FCSys.BCs.Chemical.Species.Species 'e-_an'(materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-52,-24})));
    FCSys.BCs.Chemical.Species.Species H2(materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-28,-24})));
    inner FCSys.BCs.Defaults defaults(analysis=true)
      annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    FCSys.Subregions.Reaction ORR(n_spec=4) "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{30,-10},{50,10}})));
    FCSys.BCs.Chemical.Species.Species 'e-_ca'(materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.Current,
        redeclare Modelica.Blocks.Sources.Ramp materialSpec(duration=3600e2,
          height=100*U.A)) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={14,-24})));
    FCSys.BCs.Chemical.Species.Species O2(materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={40,-24})));
    FCSys.BCs.Chemical.Species.Species H2O(materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature,
        redeclare Modelica.Blocks.Sources.Constant materialSpec(k=-2*1.20646*U.V
            /(298.15*U.K))) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={64,-24})));
    // Note:  The OCV of H2/O2 cell is 1.20646 V at 300 K, assuming the product is
    // gaseous.
  equation
    connect('e-_an'.chemical, HOR.chemical[1]) annotation (Line(
        points={{-52,-20},{-52,-10},{-40,-10},{-40,-0.666667}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(H2.chemical, HOR.chemical[3]) annotation (Line(
        points={{-28,-20},{-28,-10},{-40,-10},{-40,0.666667}},
        color={208,104,0},
        smooth=Smooth.None));
    connect('e-_ca'.chemical, ORR.chemical[1]) annotation (Line(
        points={{14,-20},{14,-10},{40,-10},{40,-0.75}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(O2.chemical, ORR.chemical[3]) annotation (Line(
        points={{40,-20},{40,0.25}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(H2O.chemical, ORR.chemical[4]) annotation (Line(
        points={{64,-20},{64,-10},{40,-10},{40,0.75}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[2], ORR.chemical[2]) annotation (Line(
        points={{-40,5.55112e-16},{8,-4.87687e-22},{8,-0.25},{40,-0.25}},
        color={208,104,0},
        smooth=Smooth.None));

    annotation (
      Diagram(graphics),
      experiment(StopTime=360000),
      experimentSetupOutput,
      Commands(file="resources/scripts/Dymola/Subregions.Examples.Cell.mos"));
  end Cell;

  model SimpleSpecies "Simple species model to test reactions"
    //extends FCSys.BaseClasses.Icons.Names.Middle;

    replaceable FCSys.Characteristics.BaseClasses.Characteristic Data
      "Characteristic data of the species"
      annotation (Dialog(group="Material properties"));
    parameter Q.Volume V=1*U.cm^3 "Volume";
    parameter Q.TemperatureAbsolute T=298.15*U.K "Temperature";
    parameter Q.PressureAbsolute p_IC=1*U.atm
      "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Scalar properties"));
    parameter Integer n_vel(
      final min=1,
      final max=3) = 1
      "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
      annotation (HideResult=true);
    parameter Boolean overrideEOS=false
      "<html>Override the equation of state with the value of &rho;<sub>IC</sub></html>"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(tab="Assumptions", compact=true),
      choices(__Dymola_checkBox=true));
    parameter Q.TemperatureAbsolute T_IC(nominal=298.15*U.K, start=defaults.T)
      "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>";
    parameter Q.AmountVolumic rho_IC(min=if overrideEOS then 0 else Modelica.Constants.small,
        start=1/Data.v_Tp(T_IC, p_IC))
      "<html>Initial volumic amount (&rho;<sub>IC</sub>)</html>";

    Q.Amount N(nominal=1*U.C, start=1*U.C) "Amount";
    Q.Pressure p(
      nominal=1*U.atm,
      start=p_IC,
      fixed=true) "Pressure";
    Q.Potential h(nominal=1*U.V) "Specific enthalpy";
    Q.Potential mu(nominal=1*U.V) "Electrochemical potential";

    FCSys.Connectors.ChemicalOutput chemical(
      final n_vel=n_vel,
      final formula=Data.formula,
      final m=Data.m) annotation (Placement(transformation(extent={{-10,-50},{
              10,-30}}), iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    outer FCSys.BCs.Defaults defaults "Default settings" annotation (Placement(
          transformation(extent={{-10,-10},{10,10}}), iconTransformation(extent
            ={{-10,90},{10,110}})));

  equation
    // Properties
    if overrideEOS then
      N = rho_IC*V;
    elseif Data.isCompressible then
      p = Data.p_Tv(T, V/N);
    else
      V = N*Data.v_Tp(T, p);
    end if;
    h = Data.h0(T);
    mu = chemical.muPerT*T;
    h = mu + T*Data.s(T, p);
    chemical.hbar = h/Data.m;
    chemical.phi = zeros(n_vel);

    // Conservation of material
    der(N)/U.s = chemical.Ndot;
    annotation (defaultComponentName="species", Icon(graphics={Rectangle(
              extent={{-100,40},{100,-40}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),Line(
              points={{-100,-40},{100,-40}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-100,-40},{-100,40},{100,40},{100,-40}},
              pattern=LinePattern.None,
              smooth=Smooth.None),Text(
              extent={{-100,-20},{100,20}},
              textString="%name",
              lineColor={0,0,0})}));
  end SimpleSpecies;

  function g_ "Gibbs potential as a function of pressure and temperature"

    extends Modelica.Icons.Function;

    input Q.PressureAbsolute p=1*U.atm "Pressure";
    input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
    input ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC
      "Choice of enthalpy reference";
    output Q.Potential g "Gibbs potential";

  protected
    function g0_i
      "Return g0 as a function of T using one of the temperature ranges, with enthalpy of formation at 25 degC"
      input Q.TemperatureAbsolute T "Temperature";
      input Integer i "index of the temperature range";
      output Q.Potential g0_i "g0";

    algorithm
      g0_i := FCSys.BaseClasses.Utilities.Poly.poly(
            T,
            {Characteristics.BaseClasses.Characteristic.b_c[i, j]*((if
          Characteristics.BaseClasses.Characteristic.specHeatCapPow + j == 0
           then ln(T) else 1/(Characteristics.BaseClasses.Characteristic.specHeatCapPow
           + j)) - (if Characteristics.BaseClasses.Characteristic.specHeatCapPow
           + j == 1 then ln(T) else 1/(Characteristics.BaseClasses.Characteristic.specHeatCapPow
           + j - 1))) - (if Characteristics.BaseClasses.Characteristic.specHeatCapPow
           + j == 1 then Characteristics.BaseClasses.Characteristic.B_c[i, 2]
           else 0) for j in 1:size(Characteristics.BaseClasses.Characteristic.b_c,
          2)},
            Characteristics.BaseClasses.Characteristic.specHeatCapPow + 1) +
        Characteristics.BaseClasses.Characteristic.B_c[i, 1]
        annotation (Inline=true);
    end g0_i;

  algorithm
    /*
    assert(T_lim_c[1] <= T and T <= T_lim_c[size(T_lim_c, 1)], "Temperature " +
    String(T/U.K) + " K is out of range for " + name + " ([" + String(T_lim_c[1]
    /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
    */
    // Note:  This is commented out so that the function can be inlined.
    // Note:  In Dymola 7.4 T_lim_c[end] can't be used instead of
    // T_lim_c[size(T_lim_c, 1)] due to:
    //    "Error, not all 'end' could be expanded."

    g := smooth(1, sum(if (Characteristics.BaseClasses.Characteristic.T_lim_c[i]
       <= T or i == 1) and (T < Characteristics.BaseClasses.Characteristic.T_lim_c[
      i + 1] or i == size(Characteristics.BaseClasses.Characteristic.T_lim_c, 1)
       - 1) then g0_i(T, i) else 0 for i in 1:size(Characteristics.BaseClasses.Characteristic.T_lim_c,
      1) - 1) + (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt0K then
      Characteristics.BaseClasses.Characteristic.Deltah0 else 0) - (if
      referenceEnthalpy == ReferenceEnthalpy.ZeroAt25degC then Characteristics.BaseClasses.Characteristic.Deltah0_f
       else 0) + Characteristics.BaseClasses.Characteristic.h_offset + sum((if
      Characteristics.BaseClasses.Characteristic.specVolPow[1] + i == 0 then ln(
      (if Characteristics.BaseClasses.Characteristic.p_min > 0 then max(p,
      Characteristics.BaseClasses.Characteristic.p_min) else p)/Characteristics.BaseClasses.Characteristic.p0)
       else (p^(Characteristics.BaseClasses.Characteristic.specVolPow[1] + i)
       - Characteristics.BaseClasses.Characteristic.p0^(Characteristics.BaseClasses.Characteristic.specVolPow[
      1] + i))/(Characteristics.BaseClasses.Characteristic.specVolPow[1] + i))*
      FCSys.BaseClasses.Utilities.Poly.poly(
        T,
        Characteristics.BaseClasses.Characteristic.b_v[i, :],
        Characteristics.BaseClasses.Characteristic.specVolPow[2] -
        Characteristics.BaseClasses.Characteristic.specVolPow[1] - i + 1) for i
       in 1:size(Characteristics.BaseClasses.Characteristic.b_v, 1)))
      annotation (
      InlineNoEvent=true,
      Inline=true,
      smoothOrder=1);
    // **Take first integral at actual pressure?

    // The first term is the integral of c_p*dT up to T with the reference
    // enthalpy at the lower bound [McBride2002, p. 2] plus T times the
    // integral of (c_p/T)*dT up to T with absolute entropy at the lower bound.
    // Both of these integrals are taken at p0.  The second polynomial is the
    // integral of v*dp from p0 to p (at T).
  end g_;

  function d2F
    "<html>Derivative of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.dF\">dF</a>()</html>"
    extends Modelica.Icons.Function;

    input Real x "Argument";
    input Real a[:] "Coefficients";
    input Integer n=0
      "Power associated with the first term (before derivatives)";
    input Real dx "Derivative of argument";
    input Real da[size(a, 1)] "Derivatives of coefficients";
    input Real d2x "Second derivative of argument";
    input Real d2a[size(a, 1)] "Second derivatives of coefficients";

    output Real d2F "Second derivative";

  algorithm
    d2F := BaseClasses.Utilities.Polynomial.df(
        x,
        a,
        n,
        dx,
        da)*dx + BaseClasses.Utilities.Polynomial.f(
        x,
        a,
        n)*d2x + BaseClasses.Utilities.Polynomial.f(
        x,
        da,
        n) + BaseClasses.Utilities.Polynomial.F(x, d2a)
      annotation (Inline=true);
  end d2F;

  function dF
    "<html>Derivative of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.F\">Y</a>()</html>"
    extends Modelica.Icons.Function;

    input Real x "Argument";
    input Real a[:] "Coefficients";
    input Integer n=0
      "Power associated with the first term (before derivative)";
    input Real dx "Derivative of argument";
    input Real da[size(a, 1)] "Derivatives of coefficients";

    output Real dF "Derivative";

  algorithm
    dF := BaseClasses.Utilities.Polynomial.f(
        x,
        a,
        n)*dx + BaseClasses.Utilities.Polynomial.F(
        x,
        da,
        n) annotation (Inline=true, derivative(order=2) = d2F);

    annotation (Documentation(info="<html>
  <p>The derivative of this function is
  <a href=\"modelica://FCSys.BaseClasses.Utilities.dpoly\">d2F</a>().</p></html>"));
  end dF;
  annotation (Commands(file="resources/scripts/units-values.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate)."));
end WorkInProgress;
