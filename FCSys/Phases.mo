within FCSys;
package Phases "Mixtures of species"
  extends Modelica.Icons.Package;

  model Gas "Gas phase"
    import Modelica.Math.BooleanVectors.countTrue;
    import FCSys;
    extends Partial(final n_spec=countTrue({inclH2,inclH2O,inclN2,inclO2}));

    // Conditionally include species.
    parameter Boolean inclH2=false "Include H2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2.Gas.Fixed H2(final n_faces) if inclH2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2 model" annotation (
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

    replaceable FCSys.Species.H2O.Gas.Fixed H2O(final n_faces) if inclH2O
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
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

    replaceable FCSys.Species.N2.Gas.Fixed N2(final n_faces) if inclN2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "N2 model" annotation (
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

    replaceable FCSys.Species.O2.Gas.Fixed O2(final n_faces) if inclO2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "O2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>O<sub>2</sub> model</html>",
        enable=inclO2),
      Placement(transformation(extent={{50,-10},{70,10}})));

    parameter Boolean inclChemical=false "Include reactions or phase change"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", compact=true),
      choices(__Dymola_checkBox=true));

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-120,-10},{-100,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-96,-34},{-76,-14}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{88,2},{108,22}}),iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{100,-10},{120,10}}), iconTransformation(extent={{70,-10},{
              90,10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{76,14},{96,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-108,-22},{-88,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{100,26},{120,46}}),iconTransformation(extent={{-60,40},{-40,
              60}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{100,-46},{120,-26}}),
          iconTransformation(extent={{70,-70},{90,-90}})));

    // Auxiliary variables (for analysis)
    output Q.PressureAbsolute p(stateSelect=StateSelect.never) = amagat.p if
      n_spec > 0 and environment.analysis "Total thermodynamic pressure";

    // Aliases
    /*
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if 
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.
  */

    Connectors.Chemical connH2(final n_trans=n_trans) if inclH2
      "Chemical connector for H2" annotation (Placement(transformation(extent={
              {-79,40},{-59,60}}), iconTransformation(extent={{-50,-50},{-30,-30}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-39,40},{-19,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));
    Connectors.Chemical connO2(final n_trans=n_trans) if inclO2
      "Chemical connector for O2" annotation (Placement(transformation(extent={
              {41,40},{61,60}}), iconTransformation(extent={{30,-50},{50,-30}})));

  protected
    Conditions.Adapters.AmagatDalton daltonAmagat if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{84,26},{64,46}})));
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              4,-48}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

  equation
    // Reactions

    // Phase change
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-28.6,5},{-28.6,50.5},{-29,50.5},{-29,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2.inter, inter) annotation (Line(
        points={{-51,-3.8},{-49,-6},{-49,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{-11,-3.8},{-9,-6},{-9,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.inter, inter) annotation (Line(
        points={{29,-3.8},{31,-6},{31,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.inter, inter) annotation (Line(
        points={{69,-3.8},{71,-6},{71,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2.direct.translational, direct.translational) annotation (Line(
        points={{-56.2,-9},{-56.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2.direct.thermal, direct.thermal) annotation (Line(
        points={{-56.2,-9},{-56.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.direct.translational, direct.translational) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.direct.thermal, direct.thermal) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.direct.translational, direct.translational) annotation (Line(
        points={{63.8,-9},{63.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.direct.thermal, direct.thermal) annotation (Line(
        points={{63.8,-9},{63.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(H2.dalton, daltonAmagat.dalton) annotation (Line(
        points={{-65,8.6},{-65,36},{70,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, daltonAmagat.dalton) annotation (Line(
        points={{-25,8.6},{-25,36},{70,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.dalton, daltonAmagat.dalton) annotation (Line(
        points={{15,8.6},{15,36},{70,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.dalton, daltonAmagat.dalton) annotation (Line(
        points={{55,8.6},{55,36},{70,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(daltonAmagat.amagat, amagat) annotation (Line(
        points={{78,36},{110,36}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2
    connect(H2.faces[facesCart[Axis.x], Side.n], xNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.x], Side.p], xPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.y], Side.n], yNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.y], Side.p], yPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{-60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.z], Side.n], zNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-48,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2.faces[facesCart[Axis.z], Side.p], zPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{-72,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-8,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-32,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // N2
    connect(N2.faces[facesCart[Axis.x], Side.n], xNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.x], Side.p], xPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.y], Side.n], yNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.y], Side.p], yPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[facesCart[Axis.z], Side.n], zNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{20,0},{32,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[facesCart[Axis.z], Side.p], zPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{8,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // O2
    connect(O2.faces[facesCart[Axis.x], Side.n], xNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.x], Side.p], xPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.y], Side.n], yNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.y], Side.p], yPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.z], Side.n], zNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{72,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.z], Side.p], zPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{48,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.chemical, connO2) annotation (Line(
        points={{51.4,5},{51.4,50.5},{51,50.5},{51,50}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(connH2, H2.chemical) annotation (Line(
        points={{-69,50},{-68.6,50},{-68.6,5}},
        color={221,23,47},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
<p>Please see the documentation of the <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-60},{
              120,60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Gas;

  model Graphite "Graphite phase"
    import assert = FCSys.Utilities.assertEval;
    import Modelica.Math.BooleanVectors.countTrue;
    import FCSys;
    extends Partial(
      final reduceTrans=false,
      final reduceThermal=true,
      final n_spec=countTrue({'inclC+','incle-'}));

    // Conditionally include species.
    parameter Boolean 'inclC+'=false "Include C+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'C+'.Graphite.Fixed 'C+'(final n_faces) if
      'inclC+' constrainedby FCSys.Species.Incompressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "C+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>C<sup>+</sup> model</html>",
        enable='inclC+'),
      Placement(transformation(extent={{-30,-10},{-10,10}})));

    parameter Boolean 'incle-'=false "Include e-" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'e-'.Graphite.Fixed 'e-'(final n_faces) if
      'incle-' constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "e- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>'e-' model</html>",
        enable='incle-'),
      Placement(transformation(extent={{10,-10},{30,10}})));

    parameter Boolean inclChemical=false "Include reactions or phase change"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", compact=true),
      choices(__Dymola_checkBox=true));

    // Reaction parameters
    // -------------------
    parameter Q.CurrentAreicAbsolute J0=1e-7*U.A/U.cm^2 if inclChemical
      "Exchange current density" annotation (Dialog(
        group="Reactions",
        enable=inclChemical,
        __Dymola_label="<html><i>J</i><sup>o</sup></html>"));
    parameter Q.CurrentAreic J_irr=0 if inclChemical
      "Irreversible current of side-reactions" annotation (Dialog(
        group="Reactions",
        enable=inclChemical,
        __Dymola_label="<html><i>J</i><sub>irr</sub></html>"));
    // These can't be outer parameters in Dymola 7.4.

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-80,-10},{-60,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-56,-34},{-36,-14}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{48,2},{68,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{60,-10},{80,10}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{36,14},{56,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-68,-22},{-48,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{60,-46},{80,-26}}), iconTransformation(
            extent={{70,-70},{90,-90}})));

    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{60,26},{80,46}}), iconTransformation(extent={{-60,40},{-40,
              60}})));

    // Aliases
    /*
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if 
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.
*/

    FCSys.Connectors.Chemical 'conne-'(final n_trans=n_trans) if 'incle-'
      "Chemical connector for e-" annotation (Placement(transformation(extent={
              {1,40},{21,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              4,-48}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Conditions.Adapters.AmagatDalton daltonAmagat if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{44,26},{24,46}})));

  equation
    // Reactions

    // Inert exchange
    connect('C+'.inter, inter) annotation (Line(
        points={{-11,-3.8},{-9,-6},{-9,-36},{70,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.inter, inter) annotation (Line(
        points={{29,-3.8},{31,-6},{31,-36},{70,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('C+'.direct.translational, direct.translational) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('C+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.translational, direct.translational) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.thermal, direct.thermal) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(daltonAmagat.amagat, amagat) annotation (Line(
        points={{38,36},{70,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(daltonAmagat.dalton, 'C+'.dalton) annotation (Line(
        points={{30,36},{-25,36},{-25,8.6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(daltonAmagat.dalton, 'e-'.dalton) annotation (Line(
        points={{30,36},{15,36},{15,8.6}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // C+
    connect('C+'.faces[facesCart[Axis.x], Side.n], xNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,5.55112e-16},{-70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[facesCart[Axis.x], Side.p], xPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{20,0},{20,5.55112e-16},{70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.y], Side.n], yNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.y], Side.p], yPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.z], Side.n], zNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-8,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[facesCart[Axis.z], Side.p], zPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-32,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // e-
    connect('e-'.faces[facesCart[Axis.x], Side.n], xNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{-20,0},{-20,5.55112e-16},{-70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.x], Side.p], xPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,5.55112e-16},{70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.y], Side.n], yNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.y], Side.p], yPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.z], Side.n], zNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{32,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('e-'.faces[facesCart[Axis.z], Side.p], zPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{8,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{11.4,5},{11.4,50.5},{11,50.5},{11,50}},
        color={221,23,47},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-60},{80,
              60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Graphite;

  model Ionomer "Ionomer phase"
    import Modelica.Math.BooleanVectors.countTrue;
    import FCSys;
    extends Partial(final n_spec=countTrue({'inclSO3-','inclH+',inclH2O}));

    parameter Q.NumberAbsolute k_common=1
      "Coupling factor for exchange among all species within the phase"
      annotation (Dialog(group="Geometry", __Dymola_label=
            "<html><i>k</i><sub>common</sub></html>"));

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

    replaceable FCSys.Species.'SO3-'.Ionomer.Fixed 'SO3-'(final n_faces) if
      'inclSO3-' constrainedby FCSys.Species.Incompressible(
      n_faces=n_faces,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "SO3- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>SO<sub>3</sub><sup>-</sup> model</html>",
        enable='inclSO3-'),
      Placement(transformation(extent={{30,-10},{50,10}})));

    parameter Boolean 'inclH+'=false "Include H+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'H+'.Ionomer.Fixed 'H+'(final n_faces) if
      'inclH+' constrainedby FCSys.Species.Incompressible(
      n_faces=n_faces,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sup>+</sup> model</html>",
        enable='inclH+'),
      Placement(transformation(extent={{-50,-10},{-30,10}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Ionomer.Fixed H2O(final n_faces, n_intra=1)
      if inclH2O constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-10,-10},{10,10}})));
    parameter Boolean inclChemical=false "Include reactions or phase change"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", compact=true),
      choices(__Dymola_checkBox=true));

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-100,-10},{-80,10}}),iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-76,-34},{-56,-14}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{68,2},{88,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{80,-10},{100,10}}),iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{56,14},{76,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-88,-22},{-68,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{80,26},{100,46}}), iconTransformation(extent={{-60,40},{-40,
              60}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{80,-40},{100,-20}}),
          iconTransformation(extent={{70,-70},{90,-90}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-19,40},{1,60}}), iconTransformation(extent={{10,-50},{30,-30}})));
    Connectors.Chemical 'connH+'(final n_trans=n_trans) if 'inclH+'
      "Chemical connector for H+" annotation (Placement(transformation(extent={
              {-59,40},{-39,60}}), iconTransformation(extent={{-30,-50},{-10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton daltonAmagat if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{64,26},{44,46}})));
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              24,-54}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Connectors.InertNode intra "Connection node for exchange within the phase"
      annotation (Placement(transformation(extent={{18,-52},{38,-32}}),
          iconTransformation(extent={{48,-76},{68,-56}})));

  equation
    // Phase change
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-8.6,5},{-8.6,50},{-9,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Reactions

    // Inert exchange
    connect('SO3-'.inter, inter) annotation (Line(
        points={{49,-3.8},{52,-6},{52,-6},{52,-6},{52,-30},{90,-30}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.translational, direct.translational) annotation (Line(
        points={{43.8,-9},{43.8,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.thermal, direct.thermal) annotation (Line(
        points={{43.8,-9},{43.8,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.inter, inter) annotation (Line(
        points={{-31,-3.8},{-28,-6},{-28,-30},{90,-30}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.translational, direct.translational) annotation (Line(
        points={{-36.2,-9},{-36.2,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-36.2,-9},{-36.2,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{9,-3.8},{12,-6},{12,-30},{90,-30}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{3.8,-9},{3.8,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{3.8,-9},{3.8,-54},{24,-54}},
        color={38,196,52},
        smooth=Smooth.None));

    connect('SO3-'.intra[1], intra.exchange) annotation (Line(
        points={{47,-7},{47,-42},{28,-42}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.intra[1], intra.exchange) annotation (Line(
        points={{-33,-7},{-33,-42},{28,-42}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.intra[1], intra.exchange) annotation (Line(
        points={{7,-7},{7,-42},{28,-42}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect('SO3-'.dalton, daltonAmagat.dalton) annotation (Line(
        points={{35,8.6},{35,36},{50,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.dalton, daltonAmagat.dalton) annotation (Line(
        points={{-45,8.6},{-45,36},{50,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, daltonAmagat.dalton) annotation (Line(
        points={{-5,8.6},{-5,36},{50,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(daltonAmagat.amagat, amagat) annotation (Line(
        points={{58,36},{90,36}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // SO3-
    connect('SO3-'.faces[facesCart[Axis.x], Side.n], xNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{-90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.x], Side.p], xPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{30,5.55112e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.y], Side.n], yNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{40,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.y], Side.p], yPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{40,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.z], Side.n], zNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{52,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('SO3-'.faces[facesCart[Axis.z], Side.p], zPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{28,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // 'H+'
    connect('H+'.faces[facesCart[Axis.x], Side.n], xNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-64,0},{-90,0},{-90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.x], Side.p], xPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-20,5.55112e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.y], Side.n], yNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-40,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.y], Side.p], yPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-40,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.z], Side.n], zNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-28,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.faces[facesCart[Axis.z], Side.p], zPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-52,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-28,6.10623e-16},{-28,5.55112e-16},{
            -90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{0,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{6.10623e-16,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{12,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-12,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{-48.6,5},{-48.6,50.5},{-49,50.5},{-49,50}},
        color={221,23,47},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html><p>Assumptions:<ol>
    <li>The water in the ionomer does not participate in the reaction (only the water vapor does).</li>
    </ol</p>
    
    <p>See <a href=\"modelica://FCSys.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> 
    for additional assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-60},{
              100,60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Ionomer;

  model Liquid "Liquid phase"
    extends Partial(final n_spec=if inclH2O then 1 else 0);

    // Conditionally include species.
    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Liquid.Fixed H2O(final n_faces) if inclH2O
      constrainedby FCSys.Species.Incompressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-20,-40},{0,-20}})));

    parameter Boolean inclChemical=false "Include reactions or phase change"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", compact=true),
      choices(__Dymola_checkBox=true));

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-60,-40},{-40,-20}}),iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-20,-80},{0,-60}}), iconTransformation(extent={{-10,-94},{
              10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{0,-20},{20,0}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{20,-40},{40,-20}}),iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{-20,0},{0,20}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-40,-60},{-20,-40}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{20,-60},{40,-40}}), iconTransformation(
            extent={{70,-70},{90,-90}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{20,20},{40,40}}), iconTransformation(extent={{-60,40},{-40,
              60}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-29,40},{-9,60}}),iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton daltonAmagat if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{20,20},{0,40}})));

  equation
    // Phase change
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-18.6,-25},{-18.6,50.5},{-19,50.5},{-19,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2O.inter, inter) annotation (Line(
        points={{-1,-33.8},{1,-36},{1,-50},{30,-50}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(daltonAmagat.amagat, amagat) annotation (Line(
        points={{14,30},{30,30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(daltonAmagat.dalton, H2O.dalton) annotation (Line(
        points={{6,30},{-15,30},{-15,-21.4}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2O
    connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{-10,-30},{-50,-30}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{-10,-30},{30,-30}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{-10,-30},{-10,-70}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{-10,-30},{-10,10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{-10,-30},{10,-10}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{-10,-30},{-30,-50}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>
    <p>Please see the documentation of the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-60,-80},{40,
              60}}), graphics));
  end Liquid;

protected
  partial model Partial "Base model for a phase"
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Middle;

    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    inner parameter Integer n_inter=0
      "Number of independent couplings with other phases"
      annotation (Dialog(connectorSizing=true),HideResult=n_spec == 0);

    // Geometric parameters
    inner parameter Q.NumberAbsolute k[Axis](
      each min=Modelica.Constants.small,
      each final nominal=1) = {1,1,1} if n_spec > 0
      "Scaling factor for diffusive transport" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><b><i>k</i></b></html>"));
    inner parameter Q.NumberAbsolute k_inter[n_inter]=ones(n_inter) if n_spec
       > 0 "Coupling factors with other phases" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>k</i><sub>inter</sub></html>"));
    parameter Integer n_faces=1 "Number of pairs of faces"
      annotation (HideResult=true);
    // This cannot be an inner/outer parameter in Dymola 2014.

    // Assumptions
    inner parameter Boolean reduceTrans=false "Same velocity for all species"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", enable=n_spec > 1),
      choices(__Dymola_checkBox=true));
    inner parameter Boolean reduceThermal=false
      "Same temperature for all species" annotation (
      HideResult=true,
      Dialog(tab="Assumptions", enable=n_spec > 1),
      choices(__Dymola_checkBox=true));

    inner Q.Volume V if n_spec > 0 "Volume";
    final inner Q.Number epsilon=V/product(L) if n_spec > 0
      "Volumetric filling ratio";

  protected
    outer parameter Q.Length L[Axis] if n_spec > 0 "Length" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Area A[Axis] if n_spec > 0 "Cross-sectional area"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    final inner Q.Length Lprime[Axis]=k*V ./ L .^ 2 if n_spec > 0
      "Effective area divided by transport length";
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer facesCart[:]
      "Face-pair indices of the Cartesian axes" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    outer parameter Boolean inclTrans[Axis]
      "true, if each component of translational momentum is included"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Boolean inclFaces[Axis]
      "true, if each pairs of faces is included" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

    outer Conditions.Environment environment "Environmental conditions";
    // This component is conditional to prevent a mathematical singularity
    // when two or more empty phases (without any species included) are
    // connected.

    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="phase",
      Documentation(info="<html><p>The scaling factor for diffusive transport (<b><i>k</i></b>) is a vector which directly affects 
    the rates of diffusion of material, transverse translational momentum, and energy of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic. Using the
    Bruggeman correction factor [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the scaling factor for diffusive transport (<b><i>k</i></b>) within a phase should be set to &epsilon;<sup>3/2</sup>
    along each axis, where &epsilon; is the volumetric filling ratio, or the ratio of the volume of the phase to the total volume of the subregion.</p>


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
  end Partial;

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
