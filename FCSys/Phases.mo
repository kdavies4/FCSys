within FCSys;
package Phases "Mixtures of species"
  extends Modelica.Icons.Package;

  model Gas "Gas phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends PartialPhase(final n_spec=countTrue({inclH2,inclH2O,inclN2,inclO2}));

    // Conditionally include species.
    parameter Boolean inclH2=false "Include H2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
        __Dymola_joinNext=true));
    replaceable FCSys.Species.H2.Gas.Fixed H2(
      final n_trans,
      final n_inter,
      final kL,
      final n_chem,
      final k_intra_Phi,
      final k_intra_Q) if inclH2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q}) "H2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub> model</html>",
        enable=inclH2),
      Placement(transformation(extent={{-70,-10},{-50,10}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Gas.Fixed H2O(
      final n_trans,
      final n_inter,
      final kL,
      final n_chem,
      final k_intra_Phi,
      final k_intra_Q) if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q}) "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-30,-10},{-10,10}})));

    parameter Boolean inclN2=false "Include N2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Nitrogen (N<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.N2.Gas.Fixed N2(
      final n_trans,
      final n_inter,
      final kL,
      final n_chem,
      final k_intra_Phi,
      final k_intra_Q) if inclN2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q}) "N2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>N<sub>2</sub> model</html>",
        enable=inclN2),
      Placement(transformation(extent={{10,-10},{30,10}})));

    parameter Boolean inclO2=false "Include O2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Oxygen (O<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.O2.Gas.Fixed O2(
      final n_trans,
      final n_inter,
      final kL,
      final n_chem,
      final k_intra_Phi,
      final k_intra_Q) if inclO2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q}) "O2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>O<sub>2</sub> model</html>",
        enable=inclO2),
      Placement(transformation(extent={{50,-10},{70,10}})));

    // Independence factors
    ExchangeParams common(k_Phi={2,2,2}, k_Q=0) if n_spec > 0
      "Among species in the phase"
      annotation (Dialog(group="Independence factors"));

    Connectors.BoundaryBus xNegative if inclTrans[Axis.x]
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-120,-10},{-100,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.BoundaryBus yNegative if inclTrans[Axis.y]
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-96,-34},{-76,-14}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.BoundaryBus zNegative if inclTrans[Axis.z]
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{88,2},{108,22}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.BoundaryBus xPositive if inclTrans[Axis.x]
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{100,-10},{120,10}}), iconTransformation(extent={{70,-10},{
              90,10}})));
    Connectors.BoundaryBus yPositive if inclTrans[Axis.y]
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{76,14},{96,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTrans[Axis.z]
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-108,-22},{-88,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-120,26},{-100,46}}), iconTransformation(extent={{70,-90},
              {90,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{100,-46},{120,-26}}),
          iconTransformation(extent={{-60,60},{-40,40}})));

    // Auxiliary variables (for analysis)
    output Q.PressureAbsolute p(stateSelect=StateSelect.never) = amagat.p if
      n_spec > 0 and environment.analysis "Total thermodynamic pressure";

    Connectors.Chemical connH2[1](each final n_trans=n_trans) if inclH2
      "Chemical connector for H2" annotation (Placement(transformation(extent={
              {-74,40},{-54,60}}), iconTransformation(extent={{-50,-50},{-30,-30}})));
    Connectors.Chemical connH2O[3](each final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-34,40},{-14,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));
    Connectors.Chemical connO2[1](each final n_trans=n_trans) if inclO2
      "Chemical connector for O2" annotation (Placement(transformation(extent={
              {46,40},{66,60}}), iconTransformation(extent={{30,-50},{50,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-100,26},{-80,46}})));
    Connectors.InertNode commonExch
      "Connector for exchange among species in the phase"
      annotation (Placement(transformation(extent={{76,-58},{96,-38}})));

  equation
    // Chemical exchange
    connect(O2.chemical, connO2) annotation (Line(
        points={{56,9},{56,50}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(H2.chemical, connH2) annotation (Line(
        points={{-64,9},{-64,50}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-24,9},{-24,50}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2.inter, inter) annotation (Line(
        points={{-51,-3.8},{-51,-36},{110,-36}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{-11,-3.8},{-11,-36},{110,-36}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(N2.inter, inter) annotation (Line(
        points={{29,-3.8},{29,-36},{110,-36}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(O2.inter, inter) annotation (Line(
        points={{69,-3.8},{69,-36},{110,-36}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.intra[1], commonExch.node) annotation (Line(
        points={{-56,-9},{-56,-48},{86,-48}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.intra[1], commonExch.node) annotation (Line(
        points={{-16,-9},{-16,-48},{86,-48}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(N2.intra[1], commonExch.node) annotation (Line(
        points={{24,-9},{24,-48},{86,-48}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(O2.intra[1], commonExch.node) annotation (Line(
        points={{64,-9},{64,-48},{86,-48}},
        color={221,23,47},
        smooth=Smooth.None));

    // Mixing
    connect(H2.dalton, amagatDalton.dalton) annotation (Line(
        points={{-69,4},{-69,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, amagatDalton.dalton) annotation (Line(
        points={{-29,4},{-29,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.dalton, amagatDalton.dalton) annotation (Line(
        points={{11,4},{11,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.dalton, amagatDalton.dalton) annotation (Line(
        points={{51,4},{51,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-94,36},{-110,36}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2
    connect(H2.boundaries[transCart[Axis.x], Side.n], xNegative.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.boundaries[transCart[Axis.x], Side.p], xPositive.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.boundaries[transCart[Axis.y], Side.n], yNegative.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{-60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.boundaries[transCart[Axis.y], Side.p], yPositive.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{-60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.boundaries[transCart[Axis.z], Side.n], zNegative.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{-48,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2.boundaries[transCart[Axis.z], Side.p], zPositive.H2) annotation
      (Line(
        points={{-60,6.10623e-016},{-72,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.boundaries[transCart[Axis.x], Side.n], xNegative.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.x], Side.p], xPositive.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.n], yNegative.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{-20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.p], yPositive.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{-20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.z], Side.n], zNegative.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{-8,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.boundaries[transCart[Axis.z], Side.p], zPositive.H2O)
      annotation (Line(
        points={{-20,6.10623e-016},{-32,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // N2
    connect(N2.boundaries[transCart[Axis.x], Side.n], xNegative.N2) annotation
      (Line(
        points={{20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.boundaries[transCart[Axis.x], Side.p], xPositive.N2) annotation
      (Line(
        points={{20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.boundaries[transCart[Axis.y], Side.n], yNegative.N2) annotation
      (Line(
        points={{20,6.10623e-016},{20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.boundaries[transCart[Axis.y], Side.p], yPositive.N2) annotation
      (Line(
        points={{20,6.10623e-016},{20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.boundaries[transCart[Axis.z], Side.n], zNegative.N2) annotation
      (Line(
        points={{20,6.10623e-016},{20,0},{32,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.boundaries[transCart[Axis.z], Side.p], zPositive.N2) annotation
      (Line(
        points={{20,6.10623e-016},{8,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // O2
    connect(O2.boundaries[transCart[Axis.x], Side.n], xNegative.O2) annotation
      (Line(
        points={{60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.boundaries[transCart[Axis.x], Side.p], xPositive.O2) annotation
      (Line(
        points={{60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.boundaries[transCart[Axis.y], Side.n], yNegative.O2) annotation
      (Line(
        points={{60,6.10623e-016},{60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.boundaries[transCart[Axis.y], Side.p], yPositive.O2) annotation
      (Line(
        points={{60,6.10623e-016},{60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.boundaries[transCart[Axis.z], Side.n], zNegative.O2) annotation
      (Line(
        points={{60,6.10623e-016},{72,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.boundaries[transCart[Axis.z], Side.p], zPositive.O2) annotation
      (Line(
        points={{60,6.10623e-016},{48,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>    
    <p>It is usually appropriate to assume that species
    exist at the same temperature within a phase.  The time constants for thermal exchange (&tau;<sub><i>Q</i>E</sub>)
    are usually much shorter than the time span of interest due to the very small coupling
    resistances.  If this is the case, set <code>commonExch.k_Q</code> to zero as <code>final</code>.  
    This will reduce the index of the problem, so it may be 
    necessary to eliminate some initial conditions.</p>

    <p>In bulk flow, it is usually appropriate to assume that species travel with the same velocity within
    a phase.  If this is the case, set the entries of <code>commonExch.k_Phi</code>
    to zero as <code>final</code>.  This will reduce the index of the problem, so it may be 
    necessary to eliminate some initial conditions.</p>
    
<p>Please see the documentation of the <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-60},{
              120,60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Gas;

  model Graphite "Graphite phase"
    import assert = FCSys.Utilities.assertEval;
    import Modelica.Math.BooleanVectors.countTrue;

    extends PartialPhase(final n_spec=countTrue({'inclC+','incle-'}));

    // Conditionally include species.
    parameter Boolean 'inclC+'=false "Include C+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'C+'.Graphite.Fixed 'C+'(
      final n_trans,
      final n_inter,
      final kL) if 'inclC+' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
      kL=kL,
      final k_intra_Phi={ones(n_trans)},
      final k_intra_Q={0}) "C+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>C<sup>+</sup> model</html>",
        enable='inclC+'),
      Placement(transformation(extent={{-30,-20},{-10,0}})));

    parameter Boolean 'incle-'=false "Include e-" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
        __Dymola_joinNext=true));
    replaceable FCSys.Species.'e-'.Graphite.Fixed 'e-'(final n_trans, final
        n_inter) if 'incle-' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=0,
      final k_intra_Phi={ones(n_trans)},
      final k_intra_Q={0}) "e- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>'e-' model</html>",
        enable='incle-'),
      Placement(transformation(extent={{10,-20},{30,0}})));

    parameter Boolean 'incle-Transfer'=false "Include electron transfer"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        enable='incle-',
        tab="Assumptions",
        compact=true));

    parameter Boolean inclDL=false "Include the double-layer capacitance"
      annotation (
      HideResult=true,
      Dialog(
        enable='incle-Transfer',
        tab="Assumptions",
        compact=true),
      choices(__Dymola_checkBox=true));

    Chemistry.Electrochemistry.ElectronTransfer 'e-Transfer'(redeclare
        constant Integer n_trans=1) if 'incle-' and 'incle-Transfer'
      "Electron transfer" annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=90,
          origin={16,34})));
    Chemistry.Electrochemistry.DoubleLayer doubleLayer(redeclare constant
        Integer n_trans=1) if 'incle-' and 'incle-Transfer' and inclDL
      "Electrolytic double layer" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-20,34})));
    // Note:  n_trans must be constant in Dymola 2014 to prevent errors
    // such as "Failed to expand the variable
    // subregion.graphite.'e-Transfer'.negative.phi".  The setting of
    // n_trans=1 must be manually changed at instantiation if additional transport
    // axes are enabled.

    Connectors.BoundaryBus xNegative if inclTrans[Axis.x]
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-80,-20},{-60,0}}), iconTransformation(extent={{-90,-10},{
              -70,10}})));
    Connectors.BoundaryBus yNegative if inclTrans[Axis.y]
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-56,-44},{-36,-24}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.BoundaryBus zNegative if inclTrans[Axis.z]
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{48,-8},{68,12}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.BoundaryBus xPositive if inclTrans[Axis.x]
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{60,-20},{80,0}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.BoundaryBus yPositive if inclTrans[Axis.y]
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{36,4},{56,24}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTrans[Axis.z]
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-68,-32},{-48,-12}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{60,-56},{80,-36}}), iconTransformation(
            extent={{-60,60},{-40,40}})));

    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-80,14},{-60,34}}), iconTransformation(extent={{70,-90},{
              90,-70}})));

    Connectors.Chemical 'conne-'[1](each final n_trans=n_trans) if 'incle-'
       and 'incle-Transfer' "Chemical connector for e-" annotation (Placement(
          transformation(extent={{6,40},{26,60}}), iconTransformation(extent={{
              -10,-50},{10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
    Connectors.InertNode commonExch
      "Connector for exchange among all species in the phase"
      annotation (Placement(transformation(extent={{36,-68},{56,-48}})));

  equation
    // Chemical exchange
    connect('e-Transfer'.negative, doubleLayer.negative) annotation (Line(
        points={{16,28},{-20,28}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-Transfer'.positive, doubleLayer.positive) annotation (Line(
        points={{16,40},{-20,40}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-Transfer'.negative, 'e-'.chemical[1]) annotation (Line(
        points={{16,28},{16,-1}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-Transfer'.positive, 'conne-'[1]) annotation (Line(
        points={{16,40},{16,50}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('C+'.inter, inter) annotation (Line(
        points={{-11,-13.8},{-11,-46},{70,-46}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('C+'.intra[1], commonExch.node) annotation (Line(
        points={{-16,-19},{-16,-58},{46,-58}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.intra[1], commonExch.node) annotation (Line(
        points={{24,-19},{24,-58},{46,-58}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-Transfer'.intra, 'C+'.intra[1]) annotation (Line(
        points={{12,34},{-16,34},{-16,-19}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(doubleLayer.intra, 'e-Transfer'.intra) annotation (Line(
        points={{-16,34},{12,34}},
        color={221,23,47},
        smooth=Smooth.None));

    // Mixing
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-54,24},{-70,24}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'C+'.dalton) annotation (Line(
        points={{-46,24},{-29,24},{-29,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'e-'.dalton) annotation (Line(
        points={{-46,24},{11,24},{11,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(doubleLayer.amagat, amagat) annotation (Line(
        points={{-20,34},{-62,34},{-62,24},{-70,24}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // C+
    connect('C+'.boundaries[transCart[Axis.x], Side.n], xNegative.'C+')
      annotation (Line(
        points={{-20,-10},{-70,-10}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.boundaries[transCart[Axis.x], Side.p], xPositive.'C+')
      annotation (Line(
        points={{-20,-10},{70,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.boundaries[transCart[Axis.y], Side.n], yNegative.'C+')
      annotation (Line(
        points={{-20,-10},{-20,-34},{-46,-34}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.boundaries[transCart[Axis.y], Side.p], yPositive.'C+')
      annotation (Line(
        points={{-20,-10},{-20,14},{46,14}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.boundaries[transCart[Axis.z], Side.n], zNegative.'C+')
      annotation (Line(
        points={{-20,-10},{-8,2},{58,2}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.boundaries[transCart[Axis.z], Side.p], zPositive.'C+')
      annotation (Line(
        points={{-20,-10},{-32,-22},{-58,-22}},
        color={127,127,127},
        smooth=Smooth.None));

    // e-
    connect('e-'.boundaries[transCart[Axis.x], Side.n], xNegative.'e-')
      annotation (Line(
        points={{20,-10},{-70,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.boundaries[transCart[Axis.x], Side.p], xPositive.'e-')
      annotation (Line(
        points={{20,-10},{70,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.boundaries[transCart[Axis.y], Side.n], yNegative.'e-')
      annotation (Line(
        points={{20,-10},{20,-34},{-46,-34}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.boundaries[transCart[Axis.y], Side.p], yPositive.'e-')
      annotation (Line(
        points={{20,-10},{20,14},{46,14}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.boundaries[transCart[Axis.z], Side.n], zNegative.'e-')
      annotation (Line(
        points={{20,-10},{32,2},{58,2}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('e-'.boundaries[transCart[Axis.z], Side.p], zPositive.'e-')
      annotation (Line(
        points={{20,-10},{8,-22},{-58,-22}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html><p>From a physical standpoint, it is probably best 
    considered its own phase, but it is modeled as a part of the graphite phase in order to reduce the number of connections.  It is assumed to 
    store no material, translational momentum, or thermal energy.  It only stores electrical energy due to a charge difference.</p>
    
    <p>C<sup>+</sup> and e<sup>-</sup> are assumed to have the same temperature.  
    C<sup>+</sup> should be used to initialize the temperature.</p>
    
    <p>See <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, please see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-80},{80,
              60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Graphite;

  model Ionomer "Ionomer phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends PartialPhase(final n_spec=countTrue({'inclSO3-','inclH+',inclH2O}),
        n_inter=1);

    // Conditionally include species.
    parameter Boolean 'inclSO3-'=false
      "Include C19HF37O5S- (abbreviated as SO3-)" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label=
            "<html>Nafion sulfonate minus (abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

        __Dymola_joinNext=true));

    replaceable FCSys.Species.'SO3-'.Ionomer.Fixed 'SO3-'(
      final n_trans,
      final n_inter,
      final kL,
      final k_intra_Phi,
      final k_intra_Q) if 'inclSO3-' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      n_intra=3,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans],protonsSolid.k_Phi[cartTrans],
          waterSolid.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q,protonsSolid.k_Q,waterSolid.k_Q}) "SO3- model"
      annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>SO<sub>3</sub><sup>-</sup> model</html>",
        enable='inclSO3-'),
      Placement(transformation(extent={{30,10},{50,30}})));

    parameter Boolean 'inclH+'=false "Include H+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'H+'.Ionomer.Fixed 'H+'(
      final n_trans,
      final n_inter,
      final n_chem,
      final kL,
      final k_intra_Phi,
      final k_intra_Q) if 'inclH+' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=3,
      n_inter=n_inter,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans],protonsSolid.k_Phi[cartTrans],
          protonsWater.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q,protonsSolid.k_Q,protonsWater.k_Q}) "H+ model"
      annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sup>+</sup> model</html>",
        enable='inclH+'),
      Placement(transformation(extent={{-50,10},{-30,30}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Ionomer.Fixed H2O(
      final n_trans,
      final n_inter,
      final n_chem,
      final kL,
      final k_intra_Phi,
      final k_intra_Q) if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=3,
      n_inter=0,
      kL=kL,
      k_intra_Phi={common.k_Phi[cartTrans],protonsWater.k_Phi[cartTrans],
          waterSolid.k_Phi[cartTrans]},
      k_intra_Q={common.k_Q,protonsWater.k_Q,waterSolid.k_Q}) "H2O model"
      annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-10,10},{10,30}})));

    Connectors.BoundaryBus xNegative if inclTrans[Axis.x]
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-100,10},{-80,30}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.BoundaryBus yNegative if inclTrans[Axis.y]
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-76,-14},{-56,6}}), iconTransformation(extent={{-10,-94},{
              10,-74}})));
    Connectors.BoundaryBus zNegative if inclTrans[Axis.z]
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{68,22},{88,42}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.BoundaryBus xPositive if inclTrans[Axis.x]
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{80,10},{100,30}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.BoundaryBus yPositive if inclTrans[Axis.y]
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{56,34},{76,54}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTrans[Axis.z]
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-88,-2},{-68,18}}), iconTransformation(extent={{-90,-90},{
              -70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-100,46},{-80,66}}), iconTransformation(extent={{70,-90},{
              90,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{80,-26},{100,-6}}), iconTransformation(
            extent={{-60,60},{-40,40}})));
    Connectors.Chemical connH2O[1](each final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-14,60},{6,80}}), iconTransformation(extent={{10,-50},{30,-30}})));
    Connectors.Chemical 'connH+'[1](each final n_trans=n_trans) if 'inclH+'
      "Chemical connector for H+" annotation (Placement(transformation(extent={
              {-54,60},{-34,80}}),iconTransformation(extent={{-30,-50},{-10,-30}})));

    // Independence factors
    ExchangeParams common(k_Q=0) if n_spec > 0 "Among all species in the phase"
      annotation (Dialog(group="Exchange (click to edit)"));
    ExchangeParams protonsSolid(final k_Phi, k_Q=Modelica.Constants.inf) if
      'inclSO3-' or 'inclH+'
      "<html>Between H<sup>+</sup> and SO<sub>3</sub><sup>-</sup>)</html>"
      annotation (Dialog(group="Exchange (click to edit)"));
    ExchangeParams protonsWater(k_Phi=fill(200, 3), k_Q=Modelica.Constants.inf)
      if 'inclH+' or inclH2O
      "<html>Between H<sup>+</sup> and H<sub>2</sub>O</html>"
      annotation (Dialog(group="Exchange (click to edit)"));

    ExchangeParams waterSolid(k_Q=Modelica.Constants.inf) if 'inclSO3-' or
      inclH2O
      "<html>Between H<sub>2</sub>O and SO<sub>3</sub><sup>-</sup></html>"
      annotation (Dialog(group="Exchange (click to edit)"));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-80,46},{-60,66}})));

    Connectors.InertNode commonExch
      "Connector for exchange among all species in the phase" annotation (
        Placement(transformation(extent={{56,-38},{76,-18}}),
          iconTransformation(extent={{86,-50},{106,-30}})));
    Connectors.InertNode protonsSolidExch
      "Connector for exchange between H+ and SO3-" annotation (Placement(
          transformation(extent={{56,-50},{76,-30}}), iconTransformation(extent
            ={{48,-76},{68,-56}})));
    Connectors.InertNode protonsWaterExch
      "Connector for exchange between H+ and H2O" annotation (Placement(
          transformation(extent={{56,-62},{76,-42}}), iconTransformation(extent
            ={{48,-76},{68,-56}})));
    Connectors.InertNode waterSolidExch
      "Connector for exchange between H2O and SO3-" annotation (Placement(
          transformation(extent={{56,-74},{76,-54}}), iconTransformation(extent
            ={{48,-76},{68,-56}})));

  equation
    // Chemical exchange
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-4,29},{-4,70}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{-44,29},{-44,70}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('SO3-'.inter, inter) annotation (Line(
        points={{49,16.2},{49,-16},{90,-16}},
        color={221,23,47},
        smooth=Smooth.None));

    connect('SO3-'.intra[1], commonExch.node) annotation (Line(
        points={{44,11},{44,-28},{66,-28}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('SO3-'.intra[2], protonsSolidExch.node) annotation (Line(
        points={{44,11},{44,-40},{66,-40}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('SO3-'.intra[3], waterSolidExch.node) annotation (Line(
        points={{44,11},{44,-64},{66,-64}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.intra[1], commonExch.node) annotation (Line(
        points={{-36,11},{-36,-28},{66,-28}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.intra[2], protonsSolidExch.node) annotation (Line(
        points={{-36,11},{-36,-40},{66,-40}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.intra[3], protonsWaterExch.node) annotation (Line(
        points={{-36,11},{-36,-52},{66,-52}},
        color={221,23,47},
        smooth=Smooth.None));

    connect(H2O.intra[1], commonExch.node) annotation (Line(
        points={{4,11},{4,-28},{66,-28}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.intra[2], waterSolidExch.node) annotation (Line(
        points={{4,11},{4,-64},{66,-64}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.intra[3], protonsWaterExch.node) annotation (Line(
        points={{4,11},{4,-52},{66,-52}},
        color={221,23,47},
        smooth=Smooth.None));

    // Mixing
    connect('SO3-'.dalton, amagatDalton.dalton) annotation (Line(
        points={{31,24},{31,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.dalton, amagatDalton.dalton) annotation (Line(
        points={{-49,24},{-49,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, amagatDalton.dalton) annotation (Line(
        points={{-9,24},{-9,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-74,56},{-90,56}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // SO3-
    connect('SO3-'.boundaries[transCart[Axis.x], Side.n], xNegative.'SO3-')
      annotation (Line(
        points={{40,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.boundaries[transCart[Axis.x], Side.p], xPositive.'SO3-')
      annotation (Line(
        points={{40,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.boundaries[transCart[Axis.y], Side.n], yNegative.'SO3-')
      annotation (Line(
        points={{40,20},{40,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.boundaries[transCart[Axis.y], Side.p], yPositive.'SO3-')
      annotation (Line(
        points={{40,20},{40,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.boundaries[transCart[Axis.z], Side.n], zNegative.'SO3-')
      annotation (Line(
        points={{40,20},{52,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('SO3-'.boundaries[transCart[Axis.z], Side.p], zPositive.'SO3-')
      annotation (Line(
        points={{40,20},{28,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));
    // 'H+'
    connect('H+'.boundaries[transCart[Axis.x], Side.n], xNegative.'H+')
      annotation (Line(
        points={{-40,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.boundaries[transCart[Axis.x], Side.p], xPositive.'H+')
      annotation (Line(
        points={{-40,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.boundaries[transCart[Axis.y], Side.n], yNegative.'H+')
      annotation (Line(
        points={{-40,20},{-40,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.boundaries[transCart[Axis.y], Side.p], yPositive.'H+')
      annotation (Line(
        points={{-40,20},{-40,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.boundaries[transCart[Axis.z], Side.n], zNegative.'H+')
      annotation (Line(
        points={{-40,20},{-28,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.boundaries[transCart[Axis.z], Side.p], zPositive.'H+')
      annotation (Line(
        points={{-40,20},{-52,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.boundaries[transCart[Axis.x], Side.n], xNegative.H2O)
      annotation (Line(
        points={{0,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.x], Side.p], xPositive.H2O)
      annotation (Line(
        points={{0,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.n], yNegative.H2O)
      annotation (Line(
        points={{0,20},{0,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.p], yPositive.H2O)
      annotation (Line(
        points={{0,20},{0,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.z], Side.n], zNegative.H2O)
      annotation (Line(
        points={{0,20},{12,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.boundaries[transCart[Axis.z], Side.p], zPositive.H2O)
      annotation (Line(
        points={{0,20},{0,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html><p>Assumptions:<ol>
    <li>The water in the ionomer does not participate in the reaction (only the water vapor does).</li>
    </ol</p>
        
    <p>See <a href=\"modelica://FCSys.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> 
    for additional assumptions.
    For more information, please see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-80},{
              100,80}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Ionomer;

  model Liquid "Liquid phase"
    extends PartialPhase(final n_spec=if inclH2O then 1 else 0);

    // Conditionally include species.
    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    // Auxiliary variables (for analysis)
    output Q.PressureAbsolute p(stateSelect=StateSelect.never) = amagat.p if
      n_spec > 0 and environment.analysis "Total thermodynamic pressure";

    replaceable FCSys.Species.H2O.Liquid.Fixed H2O(
      final n_trans,
      final n_inter,
      final kL,
      final n_chem) if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
      kL=kL) "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-20,-20},{0,0}})));

    Connectors.BoundaryBus xNegative if inclTrans[Axis.x]
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-60,-20},{-40,0}}), iconTransformation(extent={{-90,-10},{
              -70,10}})));
    Connectors.BoundaryBus yNegative if inclTrans[Axis.y]
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,-60},{0,-40}}), iconTransformation(extent={{-10,-94},{
              10,-74}})));
    Connectors.BoundaryBus zNegative if inclTrans[Axis.z]
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{0,0},{20,20}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.BoundaryBus xPositive if inclTrans[Axis.x]
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.BoundaryBus yPositive if inclTrans[Axis.y]
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,20},{0,40}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTrans[Axis.z]
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-40,-40},{-20,-20}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{20,-40},{40,-20}}), iconTransformation(
            extent={{-60,60},{-40,40}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-60,0},{-40,20}}), iconTransformation(extent={{70,-90},{90,
              -70}})));
    Connectors.Chemical connH2O[1](each final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-40,40},{-20,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

    Chemistry.SurfaceTensionContact surfaceTension(
      final n_trans=n_trans,
      final V_w=V,
      final V=product(L)) if inclH2O
      "Effect of surface tension on chemical potential" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-30,30})));

  equation
    // Chemical exchange

    // Inert exchange
    connect(H2O.inter, inter) annotation (Line(
        points={{-1,-13.8},{-1,-30},{30,-30}},
        color={221,23,47},
        smooth=Smooth.None));

    // Mixing
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-34,10},{-50,10}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, H2O.dalton) annotation (Line(
        points={{-26,10},{-19,10},{-19,-6}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2O
    connect(H2O.boundaries[transCart[Axis.x], Side.n], xNegative.H2O)
      annotation (Line(
        points={{-10,-10},{-50,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.x], Side.p], xPositive.H2O)
      annotation (Line(
        points={{-10,-10},{30,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.n], yNegative.H2O)
      annotation (Line(
        points={{-10,-10},{-10,-50}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.y], Side.p], yPositive.H2O)
      annotation (Line(
        points={{-10,-10},{-10,30}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.boundaries[transCart[Axis.z], Side.n], zNegative.H2O)
      annotation (Line(
        points={{-10,-10},{10,10}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.boundaries[transCart[Axis.z], Side.p], zPositive.H2O)
      annotation (Line(
        points={{-10,-10},{-30,-30}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(surfaceTension.nonwetting, connH2O[1]) annotation (Line(
        points={{-30,36},{-30,50}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(surfaceTension.wetting, H2O.chemical[1]) annotation (Line(
        points={{-30,30},{-30,20},{-14,20},{-14,-1}},
        color={255,195,38},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
    <p>Please see the documentation of the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-60,-60},{40,
              60}}), graphics));
  end Liquid;

protected
  partial model PartialPhase "Base model for a phase"
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Middle;

    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    parameter Integer n_inter=0
      "Number of exchange connections with other phases"
      annotation (Dialog(connectorSizing=true),HideResult=true);

    // Geometric parameters
    parameter Q.NumberAbsolute k[Axis](
      each min=Modelica.Constants.small,
      each final nominal=1) = {1,1,1} if n_spec > 0
      "Length factors for diffusive transport" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><b><i>k</i></b></html>"));
    parameter Integer n_trans=1 "Number of transport axes"
      annotation (HideResult=true);
    // This cannot be an inner/outer parameter in Dymola 2014.
    inner Q.Volume V if n_spec > 0 "Volume of the phase";

    // Independence factors
    inner parameter Q.NumberAbsolute k_inter_Phi[n_inter, n_trans]=ones(n_inter,
        n_trans) if n_spec > 0 "For translational exchange with other phases"
      annotation (Dialog(group="Independence factors", __Dymola_label=
            "<html><i>k</i><sub>inter &Phi;</sub></html>"));
    inner parameter Q.NumberAbsolute k_inter_Q[n_inter]=ones(n_inter) if n_spec
       > 0 "For thermal exchange with other phases" annotation (Dialog(group=
            "Independence factors", __Dymola_label=
            "<html><i>k</i><sub>inter <i>Q</i></sub></html>"));

    // Auxiliary variables (for analysis)
    output Q.Number epsilon=V/product(L) if n_spec > 0 and environment.analysis
      "Volumetric fill fraction";

  protected
    outer parameter Q.Length L[Axis] if n_spec > 0 "Length of the subregion"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    final parameter Q.Length kL[:]=k[cartTrans] .* L[cartTrans] if n_spec > 0
      "Effective transport lengths";
    final inner Q.Area Aprime[n_trans]=fill(V, n_trans) ./ L[cartTrans] if
      n_spec > 0 "Effective cross-sectional areas";
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer transCart[:]
      "Boundary-pair indices of the Cartesian axes" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Boolean inclTrans[Axis]
      "true, if each pair of boundaries is included" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

    outer Conditions.Environment environment "Environmental conditions";

    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="phase",
      Documentation(info="<html><p>The scaling factor for diffusive transport (<b><i>k</i></b>) is a vector which directly affects 
    the resistance associated with the transport of material, transverse translational momentum, and energy of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic.  Using
    Bruggeman correction [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the factor (<b><i>k</i></b>) within a phase should be set to &epsilon;<sup>-1/2</sup>
    along each axis, where &epsilon; is the volumetric filling ratio, or the ratio of the volume of the phase to the total volume of the subregion.
    The Bruggeman factor itself increases resistance by a &epsilon;<sup>-3/2</sup>, but a factor of &epsilon;<sup>-1</sup> is included inherently.</p>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Ellipse(
            extent={{-40,100},{40,20}},
            lineColor={127,127,127},
            startAngle=30,
            endAngle=149,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{20,-4},{100,-84}},
            lineColor={127,127,127},
            startAngle=270,
            endAngle=390,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{-100,-4},{-20,-84}},
            lineColor={127,127,127},
            startAngle=149,
            endAngle=270,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Polygon(
            points={{60,-84},{-60,-84},{-94.5,-24},{-34.5,80},{34.5,80},{94.5,-24},
                {60,-84}},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Sphere,
            smooth=Smooth.None,
            fillColor={225,225,225},
            lineColor={0,0,0}),
          Line(
            points={{-60,-84.1},{60,-84.1}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{34.5,80},{94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{-34.5,80},{-94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Text(
            extent={{-100,-20},{100,20}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(graphics));
  end PartialPhase;

public
  record ExchangeParams "Independence factors for an exchange process"
    extends Modelica.Icons.Record;

    parameter Q.NumberAbsolute k_Phi[Axis]={1,1,1} "Translational" annotation (
        Evaluate=true, Dialog(__Dymola_label=
            "<html><i>k</i><sub>&Phi;</sub></html>"));
    parameter Q.NumberAbsolute k_Q=1 "Thermal" annotation (Evaluate=true,
        Dialog(__Dymola_label="<html><i>k<sub>Q</sub></i></html>"));

  end ExchangeParams;
  annotation (Documentation(info="
<html><p>The graphite, ionomer, and
liquid phases can only be used with a compressible phase (gas).</p>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
end Phases;
